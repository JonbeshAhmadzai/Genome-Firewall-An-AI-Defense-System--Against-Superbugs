"""Genome Firewall MVP Streamlit demonstration."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

from src.genome_reader.build_features import build_feature_matrix
from src.genome_reader.fasta_qc import fasta_stats
from src.genome_reader.run_amrfinder import run_amrfinder
from src.predictor.predict import load_model, predict_target
from src.explanation.llm_explainer import explain_predictions, llm_available
from src.pipeline.router import resolve_route
from targets_config import enabled_species


ROOT = Path(__file__).resolve().parent
REPORT_DIR = ROOT / "reports" / "cohort30"
MODEL_ROOT = ROOT / "models" / "cohort30"
FEATURES_PATH = ROOT / "data" / "processed" / "cohort30" / "amrfinder_features.csv"
EVIDENCE_PATH = ROOT / "data" / "processed" / "cohort30" / "amrfinder_evidence.csv"
LABELS_PATH = ROOT / "data" / "raw" / "bvbrc" / "selected" / "selected_labels.csv"


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


@st.cache_resource
def load_cached_model(path: str) -> dict:
    """Keep model artifacts in memory so repeat uploads do not reload them."""

    return load_model(Path(path))


st.set_page_config(page_title="Genome Firewall MVP", page_icon="🧬", layout="wide")
st.markdown(
    """
    <style>
    .block-container { padding-top: 2rem; padding-bottom: 3rem; }
    [data-testid="stMetric"] {
        background: linear-gradient(135deg, #f8fbff 0%, #eef6ff 100%);
        border: 1px solid #d8e7f5;
        border-radius: 12px;
        padding: 0.75rem 1rem;
    }
    .scope-badge {
        display: inline-block; padding: 0.25rem 0.6rem; margin: 0.15rem 0.25rem 0.15rem 0;
        border-radius: 999px; background: #e8f3ff; color: #14558a; font-size: 0.85rem;
    }
    .summary-card {
        background: #f8fbff; border: 1px solid #d8e7f5; border-radius: 12px;
        padding: 0.85rem 1rem; min-height: 86px; color: #16324f;
    }
    .summary-label { font-size: 0.82rem; color: #52708d; margin-bottom: 0.3rem; }
    .summary-value { font-size: 1.65rem; font-weight: 700; color: #16324f; line-height: 1.1; }
    </style>
    """,
    unsafe_allow_html=True,
)
st.title("Genome Firewall")
st.caption("Defensive antibiotic-response decision support — research prototype")
st.warning("This is decision support, not a treatment decision. Confirm every result with standard laboratory testing.")

metrics = load_csv(REPORT_DIR / "model_metrics.csv")
benchmark = load_csv(REPORT_DIR / "model_benchmark.csv")
predictions = load_csv(REPORT_DIR / "cohort_predictions.csv")
features = load_csv(FEATURES_PATH)
evidence = load_csv(EVIDENCE_PATH)
labels = load_csv(LABELS_PATH)
available_species = (
    sorted(metrics["species"].dropna().astype(str).unique())
    if "species" in metrics and not metrics.empty
    else (
        sorted(predictions["species"].dropna().astype(str).unique())
        if "species" in predictions and not predictions.empty
        else list(enabled_species(kind="bacterium"))
    )
)


def render_decision_cards(frame: pd.DataFrame) -> None:
    """Render a compact, human-readable summary without changing the raw output."""

    if frame.empty:
        return
    st.subheader("Decision summary")
    columns = st.columns(min(3, len(frame)))
    icons = {"likely to fail": "🔴", "likely to work": "🟢", "no-call": "🟡"}
    for column, (_, row) in zip(columns, frame.iterrows()):
        decision = str(row.get("prediction", "no-call"))
        confidence = float(row.get("confidence", 0.0))
        with column:
            st.metric(
                f"{icons.get(decision, '⚪')} {row.get('antibiotic', 'antibiotic')}",
                decision,
                f"confidence {confidence:.0%}",
            )
            st.caption(f"Target gate: {row.get('target_status', 'unknown')}")


def render_summary_card(label: str, value: str, caption: str = "") -> None:
    """Render a theme-independent summary card with explicit text colors."""

    st.markdown(
        f'<div class="summary-card"><div class="summary-label">{label}</div>'
        f'<div class="summary-value">{value}</div><div class="summary-label">{caption}</div></div>',
        unsafe_allow_html=True,
    )

if predictions.empty:
    st.error("No model predictions found. Run data_fetchers/train_cohort.py first.")
    st.stop()

with st.sidebar:
    st.header("MVP scope")
    st.write(", ".join(available_species))
    st.write(f"{features['genome_id'].nunique() if 'genome_id' in features else 0} genomes")
    st.write(f"{max(len(features.columns) - 1, 0)} AMR features")
    st.write("AMRFinderPlus + calibrated logistic models")
    st.markdown('<span class="scope-badge">30-genome MVP</span><span class="scope-badge">3 antibiotics</span>', unsafe_allow_html=True)
    if llm_available():
        st.success("OpenAI explanations enabled")
    else:
        st.info("OpenAI explanations disabled (optional)")

tab_upload, tab_overview, tab_predictions, tab_evidence, tab_models = st.tabs(
    ["New FASTA", "Overview", "Predictions", "AMR evidence", "Model comparison"]
)

with tab_upload:
    st.subheader("Score a reconstructed genome")
    st.caption("Input must be one quality-checked reconstructed bacterial genome. Species identification and genome reconstruction are out of scope.")
    upload_species = st.selectbox("Species configuration", available_species)
    explain_with_llm = st.checkbox(
        "Generate an optional AI explanation",
        help="The LLM receives only structured QC, prediction, and AMR evidence. It cannot change the prediction or target gate.",
    )
    upload = st.file_uploader("Upload FASTA", type=["fa", "fna", "fasta", "gz"])
    if upload is not None and st.button("Run AMRFinderPlus and score", type="primary"):
        with tempfile.TemporaryDirectory(prefix="genome-firewall-") as temporary:
            input_path = Path(temporary) / upload.name
            output_path = Path(temporary) / "amrfinder.tsv"
            input_path.write_bytes(upload.getvalue())
            try:
                configured_drugs = sorted(metrics["antibiotic"].astype(str).unique())
                first_route = resolve_route(upload_species, configured_drugs[0], model_root=MODEL_ROOT)
                if not first_route.supported:
                    st.error(first_route.unsupported_reason())
                    st.stop()
                with st.spinner("Running FASTA QC, AMRFinderPlus, and the prediction models…"):
                    qc = fasta_stats(input_path)
                    run_amrfinder(
                        input_path,
                        output_path,
                        executable="/home/becode/miniconda3/envs/genome-firewall/bin/amrfinder",
                        organism=first_route.amrfinder_organism,
                        plus=True,
                        threads=8,
                    )
                    matrix, upload_evidence = build_feature_matrix([output_path])
                st.success("Scored successfully")
                qc_columns = st.columns(4)
                qc_columns[0].metric("Contigs", f"{qc['contigs']:,}")
                qc_columns[1].metric("Assembly size", f"{qc['total_bases']:,} bp")
                qc_columns[2].metric("N50", f"{qc['n50']:,} bp")
                qc_columns[3].metric("Ambiguous bases", f"{qc['ambiguous_fraction']:.2%}")
                all_scored = []
                for upload_drug in configured_drugs:
                    route = resolve_route(upload_species, upload_drug, model_root=MODEL_ROOT)
                    if not route.supported:
                        st.warning(route.unsupported_reason())
                        continue
                    if not route.model_path.exists():
                        st.warning(f"No trained model artifact found for {upload_species} / {upload_drug}.")
                        continue
                    artifact = load_cached_model(str(route.model_path))
                    for column in artifact["feature_columns"]:
                        if column not in matrix.columns:
                            matrix[column] = 0
                    scored = predict_target(
                        artifact,
                        matrix,
                        target_present=None,
                        target_names=route.molecular_targets,
                    )
                    scored.insert(1, "antibiotic", upload_drug)
                    scored["evidence_category"] = (
                        "known resistance gene or DNA change detected"
                        if not upload_evidence.empty
                        else "no known resistance signal found"
                    )
                    all_scored.append(scored)
                scored_frame = pd.concat(all_scored, ignore_index=True)
                render_decision_cards(scored_frame)
                with st.expander("Prediction details", expanded=True):
                    st.dataframe(
                        scored_frame,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "probability_resistant": st.column_config.NumberColumn("Resistance probability", format="%.1%"),
                            "confidence": st.column_config.NumberColumn("Confidence", format="%.1%"),
                        },
                    )
                report_json = json.dumps(
                    {"qc": qc, "predictions": json.loads(scored_frame.to_json(orient="records"))},
                    indent=2,
                    default=str,
                )
                st.download_button(
                    "Download prediction report (JSON)",
                    report_json,
                    file_name="genome_firewall_prediction_report.json",
                    mime="application/json",
                )
                st.info("Target annotation was not supplied for this upload, so likely-to-work is conservatively blocked as no-call.")
                if not upload_evidence.empty:
                    st.subheader("Detected AMR evidence")
                    st.dataframe(upload_evidence, use_container_width=True, hide_index=True)
                if explain_with_llm:
                    st.subheader("AI explanation")
                    try:
                        st.write(
                            explain_predictions(
                                scored_frame,
                                species=upload_species,
                                antibiotic="all supported antibiotics",
                                qc=qc,
                                evidence=upload_evidence,
                            )
                        )
                    except RuntimeError as error:
                        st.warning(str(error))
                        st.caption("The deterministic model output above remains valid; the LLM is optional.")
            except Exception as error:
                st.error(f"Could not score this FASTA: {error}")

with tab_overview:
    genome_count = int(features["genome_id"].nunique()) if "genome_id" in features else int(predictions["genome_id"].nunique())
    drug_count = int(metrics["antibiotic"].nunique()) if "antibiotic" in metrics else int(predictions["antibiotic"].nunique())
    feature_count = max(len(features.columns) - 1, 0) if not features.empty else 0
    first, second, third, fourth = st.columns(4)
    with first:
        render_summary_card("Genomes", f"{genome_count:,}", "active cohort")
    with second:
        render_summary_card("Drug models", f"{drug_count:,}", "supported antibiotics")
    with third:
        render_summary_card("AMR features", f"{feature_count:,}", "determinant + engineered")
    with fourth:
        render_summary_card("Evidence rows", f"{len(evidence):,}", "AMRFinderPlus findings")
    scope_species = " ".join(f"<span class='scope-badge'>{species}</span>" for species in available_species)
    scope_drugs = " ".join(
        f"<span class='scope-badge'>{drug}</span>"
        for drug in sorted(predictions["antibiotic"].astype(str).unique())
    )
    st.markdown(f"**Configured model scope:** {scope_species} {scope_drugs}", unsafe_allow_html=True)
    st.subheader("Held-out metrics")
    st.dataframe(metrics, use_container_width=True, hide_index=True)
    st.subheader("Decision policy")
    st.markdown("`likely to fail` / `likely to work` require confidence ≥ 0.70. Otherwise the system returns `no-call`. Target incompatibility also blocks a likely-to-work call.")

with tab_predictions:
    antibiotic = st.selectbox("Antibiotic", sorted(predictions["antibiotic"].unique()))
    view = predictions[predictions["antibiotic"].eq(antibiotic)].copy()
    st.dataframe(
        view[["genome_id", "probability_resistant", "confidence", "prediction", "target_status", "target_gate_reason"]],
        use_container_width=True,
        hide_index=True,
    )
    st.download_button("Download prediction table", view.to_csv(index=False), f"{antibiotic}_predictions.csv", "text/csv")

with tab_evidence:
    if evidence.empty:
        st.info("No AMRFinder evidence table found.")
    else:
        genome = st.selectbox("Genome", sorted(evidence["genome_id"].astype(str).unique()))
        selected = evidence[evidence["genome_id"].astype(str).eq(genome)]
        st.dataframe(selected, use_container_width=True, hide_index=True)
        st.caption("Evidence is presented as detected AMR determinants; statistical model association is not biological causation.")

with tab_models:
    if benchmark.empty:
        st.info("Run the exploratory benchmark to compare models.")
    else:
        st.dataframe(benchmark, use_container_width=True, hide_index=True)
        chart_data = benchmark.pivot_table(index="model", columns="antibiotic", values="balanced_accuracy")
        if not chart_data.empty:
            st.subheader("Balanced accuracy comparison")
            st.bar_chart(chart_data)
        st.caption("The 30-genome MVP is too small for reliable model ranking; this table is for plumbing and demo comparison only.")

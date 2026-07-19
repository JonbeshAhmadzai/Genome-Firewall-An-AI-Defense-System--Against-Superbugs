"""Genome Firewall MVP Streamlit demonstration."""

from __future__ import annotations

from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

from src.genome_reader.build_features import build_feature_matrix
from src.genome_reader.fasta_qc import fasta_stats
from src.genome_reader.run_amrfinder import run_amrfinder
from src.predictor.predict import load_model, predict_target
from targets_config import MOLECULAR_TARGETS
from src.explanation.llm_explainer import explain_predictions, llm_available


ROOT = Path(__file__).resolve().parent
REPORT_DIR = ROOT / "reports" / "cohort30"
FEATURES_PATH = ROOT / "data" / "processed" / "cohort30" / "amrfinder_features.csv"
EVIDENCE_PATH = ROOT / "data" / "processed" / "cohort30" / "amrfinder_evidence.csv"
LABELS_PATH = ROOT / "data" / "raw" / "bvbrc" / "selected" / "selected_labels.csv"


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


st.set_page_config(page_title="Genome Firewall MVP", page_icon="🧬", layout="wide")
st.title("Genome Firewall")
st.caption("Defensive antibiotic-response decision support — research prototype")
st.warning("This is decision support, not a treatment decision. Confirm every result with standard laboratory testing.")

metrics = load_csv(REPORT_DIR / "model_metrics.csv")
benchmark = load_csv(REPORT_DIR / "model_benchmark.csv")
predictions = load_csv(REPORT_DIR / "cohort_predictions.csv")
features = load_csv(FEATURES_PATH)
evidence = load_csv(EVIDENCE_PATH)
labels = load_csv(LABELS_PATH)

if predictions.empty:
    st.error("No model predictions found. Run data_fetchers/train_cohort.py first.")
    st.stop()

with st.sidebar:
    st.header("MVP scope")
    st.write("Escherichia coli")
    st.write(f"{features['genome_id'].nunique() if 'genome_id' in features else 0} genomes")
    st.write(f"{max(len(features.columns) - 1, 0)} AMR features")
    st.write("AMRFinderPlus + calibrated logistic models")
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
                qc = fasta_stats(input_path)
                run_amrfinder(
                    input_path,
                    output_path,
                    executable="/home/becode/miniconda3/envs/genome-firewall/bin/amrfinder",
                    organism="Escherichia",
                    plus=True,
                    threads=8,
                )
                matrix, upload_evidence = build_feature_matrix([output_path])
                st.success("Scored successfully")
                st.caption(f"FASTA QC: {qc['contigs']} contigs, {qc['total_bases']:,} bases, N50 {qc['n50']:,}, ambiguous fraction {qc['ambiguous_fraction']:.4f}")
                all_scored = []
                for upload_drug in sorted(metrics["antibiotic"].unique()):
                    model_path = ROOT / "models" / "cohort30" / f"escherichia_coli__{upload_drug}.joblib"
                    artifact = load_model(model_path)
                    for column in artifact["feature_columns"]:
                        if column not in matrix.columns:
                            matrix[column] = 0
                    scored = predict_target(
                        artifact,
                        matrix,
                        target_present=None,
                        target_names=MOLECULAR_TARGETS.get(upload_drug, ()),
                    )
                    scored.insert(1, "antibiotic", upload_drug)
                    scored["evidence_category"] = (
                        "known resistance gene or DNA change detected"
                        if not upload_evidence.empty
                        else "no known resistance signal found"
                    )
                    all_scored.append(scored)
                st.dataframe(pd.concat(all_scored, ignore_index=True), use_container_width=True, hide_index=True)
                st.info("Target annotation was not supplied for this upload, so likely-to-work is conservatively blocked as no-call.")
                if not upload_evidence.empty:
                    st.subheader("Detected AMR evidence")
                    st.dataframe(upload_evidence, use_container_width=True, hide_index=True)
                if explain_with_llm:
                    st.subheader("AI explanation")
                    try:
                        st.write(
                            explain_predictions(
                                pd.concat(all_scored, ignore_index=True),
                                species="Escherichia coli",
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
    first, second, third, fourth = st.columns(4)
    first.metric("Genomes", int(features["genome_id"].nunique()))
    second.metric("Drug models", int(metrics["antibiotic"].nunique()) if not metrics.empty else 0)
    third.metric("AMR features", max(len(features.columns) - 1, 0))
    fourth.metric("Evidence rows", len(evidence))
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
        st.caption("The 30-genome MVP is too small for reliable model ranking; this table is for plumbing and demo comparison only.")

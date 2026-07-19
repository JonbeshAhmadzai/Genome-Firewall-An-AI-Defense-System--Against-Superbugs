from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st
import yaml


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.features.kmer import kmer_feature_row, read_fasta_text
from genome_firewall.models.predict import predict_feature_table


MODEL_CONFIGS = {
    "AMRFinderPlus": {
        "manifest": ROOT / "models" / "amrfinder" / "manifest.json",
        "predictions": ROOT / "reports" / "metrics" / "amrfinder_cohort_predictions.csv",
        "metrics": ROOT / "reports" / "metrics" / "amrfinder_model_metrics.csv",
        "model_card": ROOT / "reports" / "amrfinder_model_card.md",
        "upload_enabled": False,
    },
    "Temporary k-mer": {
        "manifest": ROOT / "models" / "genome_firewall" / "manifest.json",
        "predictions": ROOT / "reports" / "metrics" / "cohort_predictions.csv",
        "metrics": ROOT / "reports" / "metrics" / "model_metrics.csv",
        "model_card": ROOT / "reports" / "model_card.md",
        "upload_enabled": True,
    },
}
ANTIBIOTICS_PATH = ROOT / "configs" / "antibiotics.yaml"


st.set_page_config(page_title="Genome Firewall", layout="wide")


@st.cache_data
def load_manifest(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


@st.cache_data
def load_antibiotics() -> dict:
    return yaml.safe_load(ANTIBIOTICS_PATH.read_text(encoding="utf-8"))


@st.cache_data
def load_metrics(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


@st.cache_data
def load_cohort_predictions(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


@st.cache_data
def load_model_card(path: str) -> str:
    model_card_path = Path(path)
    if not model_card_path.exists():
        return "Model card has not been generated yet."
    return model_card_path.read_text(encoding="utf-8")


def predict_uploaded_fasta(uploaded_file, manifest: dict, manifest_path: Path) -> pd.DataFrame:
    text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    sequence = read_fasta_text(text)
    if not sequence:
        raise ValueError("No DNA sequence found in uploaded FASTA.")

    feature_type = manifest.get("feature_type", "")
    if not feature_type.startswith("temporary_4mer"):
        raise ValueError("Uploaded FASTA prediction currently requires the temporary k-mer model.")

    first_model = next(iter(manifest["models"].values()))
    feature_columns = first_model["feature_columns"]
    row = {"genome_id": uploaded_file.name}
    row.update(kmer_feature_row(sequence, feature_columns, k=4))
    features = pd.DataFrame([row])
    return predict_feature_table(features, manifest_path, ROOT)


def decision_badge(decision: str) -> str:
    if decision == "likely_to_fail":
        return "Likely to fail"
    if decision == "likely_to_work":
        return "Likely to work"
    return "No-call"


def render_prediction_table(predictions: pd.DataFrame, antibiotics: dict) -> None:
    display = predictions.copy()
    display["result"] = display["decision"].map(decision_badge)
    display["confidence"] = display["confidence"].map(lambda value: f"{value:.2f}")
    display["resistant_probability"] = display["resistant_probability"].map(lambda value: f"{value:.2f}")
    display["drug_class"] = display["antibiotic"].map(
        lambda drug: antibiotics["antibiotics"].get(drug, {}).get("class", "not configured")
    )
    if "target_gate_status" not in display.columns:
        display["target_gate_status"] = display["antibiotic"].map(
            lambda drug: antibiotics["antibiotics"].get(drug, {}).get("target_presence_rule", "not configured")
        )
    if "molecular_targets" not in display.columns:
        display["molecular_targets"] = display["antibiotic"].map(
            lambda drug: ", ".join(antibiotics["antibiotics"].get(drug, {}).get("primary_targets", []))
        )
    st.dataframe(
        display[
            [
                "antibiotic",
                "drug_class",
                "result",
                "confidence",
                "resistant_probability",
                "evidence_category",
                *([column for column in ["supporting_amr_hits", "supporting_amr_classes"] if column in display.columns]),
                "target_gate_status",
                *([column for column in ["target_present", "matched_target_markers"] if column in display.columns]),
                "molecular_targets",
            ]
        ],
        hide_index=True,
        use_container_width=True,
    )


def main() -> None:
    with st.sidebar:
        selected_model = st.radio("Model source", list(MODEL_CONFIGS.keys()), index=0)
    config = MODEL_CONFIGS[selected_model]
    manifest_path = config["manifest"]
    manifest = load_manifest(str(manifest_path))
    antibiotics = load_antibiotics()

    st.title("Genome Firewall")
    st.caption("Defensive antibiotic-response prediction from reconstructed bacterial genomes.")
    st.warning(manifest.get("safety_warning", "Research prototype only. Confirm with standard laboratory testing."))

    with st.sidebar:
        st.subheader("Model")
        st.write(f"Selected: `{selected_model}`")
        st.write(f"Feature source: `{manifest.get('feature_type', 'unknown')}`")
        st.write(f"Supported species: `{antibiotics['supported_species']}`")
        st.write(f"Models: `{', '.join(manifest.get('models', {}).keys())}`")
        st.write("Thresholds")
        st.json(manifest.get("thresholds", {}))

    tab_upload, tab_cohort, tab_metrics, tab_card = st.tabs(
        ["FASTA Prediction", "Cohort Review", "Evaluation", "Model Card"]
    )

    with tab_upload:
        uploaded = st.file_uploader("Upload reconstructed E. coli FASTA", type=["fna", "fa", "fasta", "txt"])
        if not config["upload_enabled"]:
            st.info(
                "Live FASTA upload is disabled for the AMRFinderPlus model in this Windows app. "
                "Run AMRFinderPlus in WSL, rebuild features, and review predictions in the Cohort Review tab. "
                "Use the temporary k-mer model only for upload-flow demonstration."
            )
        elif uploaded is not None:
            try:
                predictions = predict_uploaded_fasta(uploaded, manifest, manifest_path)
                render_prediction_table(predictions, antibiotics)
            except Exception as exc:
                st.error(str(exc))
        else:
            st.info("Upload a FASTA file to generate antibiotic-response predictions.")

    with tab_cohort:
        cohort = load_cohort_predictions(str(config["predictions"]))
        genome_ids = sorted(cohort["genome_id"].astype(str).unique())
        selected = st.selectbox("Genome", genome_ids)
        render_prediction_table(cohort[cohort["genome_id"].astype(str) == selected], antibiotics)

    with tab_metrics:
        metrics = load_metrics(str(config["metrics"]))
        st.dataframe(metrics, hide_index=True, use_container_width=True)
        st.write("Reliability and called-decision plots are written to `reports/figures/`.")

    with tab_card:
        st.markdown(load_model_card(str(config["model_card"])))


if __name__ == "__main__":
    main()

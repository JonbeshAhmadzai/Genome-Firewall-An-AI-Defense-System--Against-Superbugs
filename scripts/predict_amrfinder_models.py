from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.models.predict import predict_feature_table
from genome_firewall.models.target_gate import apply_target_gate


FEATURES_PATH = ROOT / "data" / "interim" / "ecoli_amrfinder_features.csv"
EVIDENCE_PATH = ROOT / "data" / "interim" / "ecoli_amrfinder_evidence.csv"
METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
TARGET_PRESENCE_PATH = ROOT / "data" / "interim" / "target_presence.csv"
MANIFEST_PATH = ROOT / "models" / "amrfinder" / "manifest.json"
OUT_PATH = ROOT / "reports" / "metrics" / "amrfinder_cohort_predictions.csv"
ANTIBIOTICS_PATH = ROOT / "configs" / "antibiotics.yaml"


def summarize_evidence(evidence: pd.DataFrame, antibiotics: dict) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame(
            columns=[
                "genome_id",
                "antibiotic",
                "supporting_amr_hits",
                "supporting_amr_classes",
                "all_amr_hits",
            ]
        )
    evidence = evidence.copy()
    evidence["genome_id"] = evidence["genome_id"].astype(str)
    all_hits = evidence.groupby("genome_id").agg(
        all_amr_hits=(
            "element_symbol",
            lambda values: ", ".join(sorted({str(value) for value in values if str(value) != "nan"})[:12]),
        )
    ).reset_index()

    rows = []
    for antibiotic, metadata in antibiotics.get("antibiotics", {}).items():
        relevant_classes = {str(value).upper() for value in metadata.get("relevant_amr_classes", [])}
        relevant = evidence[evidence["class"].astype(str).str.upper().isin(relevant_classes)]
        if relevant.empty:
            continue
        grouped = relevant.groupby("genome_id").agg(
            supporting_amr_hits=(
                "element_symbol",
                lambda values: ", ".join(sorted({str(value) for value in values if str(value) != "nan"})[:12]),
            ),
            supporting_amr_classes=(
                "class",
                lambda values: ", ".join(
                    sorted({str(value) for value in values if str(value) and str(value) != "nan"})[:8]
                ),
            ),
        ).reset_index()
        grouped["antibiotic"] = antibiotic
        rows.append(grouped)

    if rows:
        summary = pd.concat(rows, ignore_index=True)
    else:
        summary = pd.DataFrame(columns=["genome_id", "antibiotic", "supporting_amr_hits", "supporting_amr_classes"])
    summary = summary.merge(all_hits, on="genome_id", how="left")
    return summary


def main() -> None:
    features = pd.read_csv(FEATURES_PATH)
    evidence = pd.read_csv(EVIDENCE_PATH) if EVIDENCE_PATH.exists() else pd.DataFrame()
    antibiotics = yaml.safe_load(ANTIBIOTICS_PATH.read_text(encoding="utf-8"))
    metadata = pd.read_csv(METADATA_PATH)
    metadata["genome_id"] = metadata["genome_id"].astype(str)
    features["genome_id"] = features["genome_id"].astype(str)
    features = features[features["genome_id"].isin(set(metadata["genome_id"]))].copy()
    predictions = predict_feature_table(features, MANIFEST_PATH, ROOT)
    target_presence = pd.read_csv(TARGET_PRESENCE_PATH) if TARGET_PRESENCE_PATH.exists() else None
    predictions = apply_target_gate(predictions, metadata, antibiotics, target_presence)
    predictions["genome_id"] = predictions["genome_id"].astype(str)
    evidence_summary = summarize_evidence(evidence, antibiotics)
    report = predictions.merge(evidence_summary, on=["genome_id", "antibiotic"], how="left")
    report["supporting_amr_hits"] = report["supporting_amr_hits"].fillna("")
    report["supporting_amr_classes"] = report["supporting_amr_classes"].fillna("")
    report["all_amr_hits"] = report["all_amr_hits"].fillna("")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(OUT_PATH, index=False)
    print(f"Wrote {OUT_PATH}")
    print(report.head(12).to_string(index=False))


if __name__ == "__main__":
    main()

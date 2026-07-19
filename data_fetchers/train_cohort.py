#!/usr/bin/env python3
"""Train and evaluate the calibrated per-antibiotic MVP models."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.predictor.predict import load_model, predict_target
from src.predictor.train import prepare_target_frame, train_target
from src.validation.evaluate import evaluate_artifact
from targets_config import MOLECULAR_TARGETS


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=Path("data/raw/bvbrc/selected/selected_labels.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/cohort20_genome_metadata.csv"))
    parser.add_argument("--groups", type=Path, default=Path("data/processed/cohort20/genome_groups.csv"))
    parser.add_argument("--features", type=Path, default=Path("data/processed/cohort20/amrfinder_features.csv"))
    parser.add_argument("--evidence", type=Path, default=Path("data/processed/cohort20/amrfinder_evidence.csv"))
    parser.add_argument("--target-presence", type=Path, default=Path("data/processed/cohort20_target_presence.csv"))
    parser.add_argument("--model-dir", type=Path, default=Path("models/cohort20"))
    parser.add_argument("--report-dir", type=Path, default=Path("reports/cohort20"))
    args = parser.parse_args()

    labels = pd.read_csv(args.labels, dtype=str)
    features = pd.read_csv(args.features, dtype={"genome_id": str})
    groups = pd.read_csv(args.groups, dtype=str)
    metadata = pd.read_csv(args.metadata, dtype=str)
    target_presence = pd.read_csv(args.target_presence, dtype=str)
    args.model_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    metrics: list[dict] = []
    predictions: list[pd.DataFrame] = []
    manifest: dict[str, object] = {"models": {}, "skipped": {}}

    for species, antibiotic in labels[["species", "antibiotic"]].drop_duplicates().itertuples(index=False):
        if antibiotic not in MOLECULAR_TARGETS:
            continue
        try:
            frame = prepare_target_frame(labels, features, groups, species=species, antibiotic=antibiotic)
            receipt = train_target(frame, output_dir=args.model_dir, species=species, antibiotic=antibiotic)
            artifact = load_model(Path(receipt["model_path"]))
            split_frame = pd.read_csv(receipt["split_path"])
            slug = Path(receipt["model_path"]).stem
            evaluation = evaluate_artifact(artifact, split_frame, output_dir=args.report_dir / slug)
            row = {"species": species, "antibiotic": antibiotic, **evaluation["metrics"]}
            metrics.append(row)
            target_rows = target_presence[target_presence["antibiotic"].eq(antibiotic)].set_index("genome_id")
            target_values = [
                (str(target_rows.loc[str(genome_id), "target_present"]).lower() == "true")
                if str(genome_id) in target_rows.index
                else None
                for genome_id in features["genome_id"].astype(str)
            ]
            predicted = predict_target(
                artifact,
                features,
                target_present=target_values,
                target_names=MOLECULAR_TARGETS[antibiotic],
            )
            predicted.insert(1, "species", species)
            predicted.insert(2, "antibiotic", antibiotic)
            evidence_ids = set(
                pd.read_csv(args.evidence, dtype=str)
                .query("genome_id.notna()")
                .loc[lambda frame: frame["genome_id"].astype(str).isin(features["genome_id"].astype(str)), "genome_id"]
                .astype(str)
            )
            predicted["evidence_category"] = predicted["genome_id"].astype(str).map(
                lambda genome_id: "known resistance gene or DNA change detected"
                if genome_id in evidence_ids
                else "no known resistance signal found"
            )
            predicted["safety_warning"] = "Research prototype only. Confirm with standard laboratory testing."
            predictions.append(predicted)
            manifest["models"][antibiotic] = receipt
        except Exception as error:
            manifest["skipped"][antibiotic] = str(error)
            print(f"Skipping {antibiotic}: {error}")

    pd.DataFrame(metrics).to_csv(args.report_dir / "model_metrics.csv", index=False)
    report = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()
    report.to_csv(args.report_dir / "cohort_predictions.csv", index=False)
    (args.report_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"Trained {len(metrics)} models; wrote {args.report_dir}")
    return 0 if metrics else 1


if __name__ == "__main__":
    raise SystemExit(main())

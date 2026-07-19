from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.models.training import load_model_table, train_models, write_manifest
from genome_firewall.models.decisions import DecisionThresholds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", type=Path, default=ROOT / "data" / "processed" / "ecoli_cohort_labels.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data" / "processed" / "ecoli_genome_metadata.csv")
    parser.add_argument("--features", type=Path, default=ROOT / "data" / "interim" / "ecoli_kmer_features.csv")
    parser.add_argument("--feature-prefix", default="kmer_")
    parser.add_argument("--feature-type", default="temporary_4mer_frequency")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "models" / "genome_firewall")
    parser.add_argument("--metrics-out", type=Path, default=ROOT / "reports" / "metrics" / "model_metrics.csv")
    parser.add_argument("--predictions-out", type=Path, default=ROOT / "reports" / "metrics" / "heldout_predictions.csv")
    parser.add_argument("--manifest-out", type=Path, default=ROOT / "models" / "genome_firewall" / "manifest.json")
    parser.add_argument("--split-audit-out", type=Path, default=ROOT / "reports" / "metrics" / "split_audit.csv")
    parser.add_argument("--min-class-count", type=int, default=8)
    parser.add_argument("--likely-to-work-max", type=float, default=0.30)
    parser.add_argument("--likely-to-fail-min", type=float, default=0.70)
    args = parser.parse_args()

    if args.likely_to_work_max >= args.likely_to_fail_min:
        raise SystemExit("--likely-to-work-max must be lower than --likely-to-fail-min")

    table, feature_columns = load_model_table(
        labels_path=args.labels,
        metadata_path=args.metadata,
        features_path=args.features,
        feature_prefix=args.feature_prefix,
    )
    metrics, predictions, manifest, split_audit = train_models(
        table=table,
        feature_columns=feature_columns,
        output_dir=args.output_dir,
        feature_type=args.feature_type,
        min_class_count=args.min_class_count,
        thresholds=DecisionThresholds(
            likely_to_work_max=args.likely_to_work_max,
            likely_to_fail_min=args.likely_to_fail_min,
        ),
    )

    args.metrics_out.parent.mkdir(parents=True, exist_ok=True)
    args.predictions_out.parent.mkdir(parents=True, exist_ok=True)
    args.split_audit_out.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.metrics_out, index=False)
    predictions.to_csv(args.predictions_out, index=False)
    split_audit.to_csv(args.split_audit_out, index=False)
    write_manifest(manifest, args.manifest_out, ROOT)

    print(f"Wrote {args.metrics_out}")
    print(f"Wrote {args.predictions_out}")
    print(f"Wrote {args.split_audit_out}")
    print(f"Wrote {args.manifest_out}")


if __name__ == "__main__":
    main()

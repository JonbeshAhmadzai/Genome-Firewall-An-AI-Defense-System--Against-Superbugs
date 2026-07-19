from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.models.predict import predict_feature_table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=ROOT / "data" / "interim" / "ecoli_kmer_features.csv")
    parser.add_argument("--manifest", type=Path, default=ROOT / "models" / "genome_firewall" / "manifest.json")
    parser.add_argument("--out", type=Path, default=ROOT / "reports" / "metrics" / "cohort_predictions.csv")
    args = parser.parse_args()

    features = pd.read_csv(args.features)
    predictions = predict_feature_table(features, args.manifest, ROOT)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")
    print(predictions.head(12).to_string(index=False))


if __name__ == "__main__":
    main()


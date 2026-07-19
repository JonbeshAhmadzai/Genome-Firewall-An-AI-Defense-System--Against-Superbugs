#!/usr/bin/env python3
"""Exploratory model comparison on the same grouped test split."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import joblib
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.predictor.grouped_split import grouped_three_way_split
from src.predictor.train import prepare_target_frame
from src.validation.scorecard import scorecard


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=Path("data/raw/bvbrc/selected/selected_labels.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/cohort20_genome_metadata.csv"))
    parser.add_argument("--groups", type=Path, default=Path("data/processed/cohort20/genome_groups.csv"))
    parser.add_argument("--features", type=Path, default=Path("data/processed/cohort20/amrfinder_features.csv"))
    parser.add_argument("--out", type=Path, default=Path("reports/cohort20/model_benchmark.csv"))
    parser.add_argument("--model-dir", type=Path, default=Path("models/cohort20/benchmarks"))
    parser.add_argument("--test-size", type=float, default=0.20, help="Fraction reserved for the grouped test split")
    parser.add_argument("--calibration-size", type=float, default=0.20, help="Fraction reserved for grouped probability calibration")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for the grouped split")
    args = parser.parse_args()
    labels = pd.read_csv(args.labels, dtype=str)
    metadata = pd.read_csv(args.metadata, dtype=str)
    groups = pd.read_csv(args.groups, dtype=str)
    features = pd.read_csv(args.features, dtype={"genome_id": str})
    rows: list[dict] = []
    args.model_dir.mkdir(parents=True, exist_ok=True)
    estimators = {
        "logistic_regression": LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
        "random_forest": RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42),
        "extra_trees": ExtraTreesClassifier(n_estimators=200, class_weight="balanced", random_state=42),
        "xgboost": XGBClassifier(
            n_estimators=100,
            max_depth=2,
            learning_rate=0.05,
            subsample=0.9,
            colsample_bytree=0.9,
            objective="binary:logistic",
            eval_metric="logloss",
            tree_method="hist",
            random_state=42,
            n_jobs=2,
        ),
    }
    for species, antibiotic in labels[["species", "antibiotic"]].drop_duplicates().itertuples(index=False):
        try:
            frame = prepare_target_frame(labels, features, groups, species=species, antibiotic=antibiotic)
            assignment = grouped_three_way_split(
                frame,
                test_size=args.test_size,
                calibration_size=args.calibration_size,
                random_state=args.seed,
            )
            split = assignment.frame
            excluded = {"genome_id", "species", "antibiotic", "phenotype", "y", "homology_group_id", "split", "_original_index"}
            columns = [column for column in split.columns if column not in excluded and pd.api.types.is_numeric_dtype(split[column])]
            train = split[split.split.eq("train")]
            test = split[split.split.eq("test")]
            if train.y.nunique() < 2 or test.y.nunique() < 2:
                continue
            for name, estimator in estimators.items():
                estimator.fit(train[columns], train.y.astype(int))
                probability = estimator.predict_proba(test[columns])[:, 1]
                rows.append({"antibiotic": antibiotic, "model": name, **scorecard(test.y, probability)})
                slug = f"{species.lower().replace(' ', '_')}__{antibiotic}__{name}"
                joblib.dump(
                    {
                        "model": estimator,
                        "feature_columns": columns,
                        "species": species,
                        "antibiotic": antibiotic,
                        "model_name": name,
                        "warning": "Exploratory MVP benchmark; not validated for clinical use.",
                    },
                    args.model_dir / f"{slug}.joblib",
                )
        except ValueError:
            continue
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(rows)} comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

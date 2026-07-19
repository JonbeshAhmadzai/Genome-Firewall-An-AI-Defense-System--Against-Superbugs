from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.models.decisions import DecisionThresholds
from genome_firewall.models.metrics import binary_metrics, no_call_metrics
from genome_firewall.models.training import grouped_train_test_split, load_model_table


def build_logistic(random_state: int) -> CalibratedClassifierCV:
    base = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=2000,
                    solver="liblinear",
                    random_state=random_state,
                ),
            ),
        ]
    )
    return CalibratedClassifierCV(estimator=base, method="sigmoid", cv=3)


def feature_sets(feature_columns: list[str]) -> dict[str, list[str]]:
    return {
        "all_refined": feature_columns,
        "family_only": [column for column in feature_columns if column.startswith("refined_family_")],
        "class_only": [column for column in feature_columns if column.startswith("refined_class_")],
        "subclass_only": [column for column in feature_columns if column.startswith("refined_subclass_")],
        "summary_only": [column for column in feature_columns if column.startswith("refined_summary_")],
        "family_plus_summary": [
            column
            for column in feature_columns
            if column.startswith("refined_family_") or column.startswith("refined_summary_")
        ],
        "class_plus_family": [
            column
            for column in feature_columns
            if column.startswith("refined_class_") or column.startswith("refined_family_")
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Ablate refined AMRFinder feature groups with one fixed model.")
    parser.add_argument("--labels", type=Path, default=ROOT / "data" / "processed" / "ecoli_cohort_labels.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data" / "processed" / "ecoli_genome_metadata.csv")
    parser.add_argument(
        "--features",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refined_features.csv",
    )
    parser.add_argument("--feature-prefix", default="refined_")
    parser.add_argument("--out", type=Path, default=ROOT / "reports" / "metrics" / "refined_feature_ablation_logistic.csv")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--test-size", type=float, default=0.25)
    args = parser.parse_args()

    table, feature_columns = load_model_table(args.labels, args.metadata, args.features, args.feature_prefix)
    thresholds = DecisionThresholds()
    rows: list[dict] = []

    for antibiotic, antibiotic_rows in table.groupby("antibiotic"):
        train, test = grouped_train_test_split(antibiotic_rows, args.random_state, args.test_size)
        for set_name, columns in feature_sets(feature_columns).items():
            if not columns:
                continue
            model = build_logistic(args.random_state)
            model.fit(train[columns], train["y"])
            probabilities = model.predict_proba(test[columns])[:, 1]
            metrics = binary_metrics(test["y"].to_numpy(), probabilities)
            metrics.update(no_call_metrics(test["y"].to_numpy(), probabilities, thresholds))
            rows.append(
                {
                    "antibiotic": antibiotic,
                    "feature_set": set_name,
                    "n_features": len(columns),
                    "test_rows": len(test),
                    **metrics,
                }
            )

    results = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")
    print(results.pivot(index="antibiotic", columns="feature_set", values="balanced_accuracy").round(3).to_string())
    print("\nBest feature set per antibiotic:")
    best = (
        results.sort_values(["antibiotic", "balanced_accuracy", "brier_score"], ascending=[True, False, True])
        .groupby("antibiotic")
        .head(1)
    )
    print(
        best[
            [
                "antibiotic",
                "feature_set",
                "n_features",
                "balanced_accuracy",
                "resistant_recall",
                "susceptible_recall",
                "brier_score",
                "no_call_rate",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()

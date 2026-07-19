from __future__ import annotations

import argparse
import sys
from pathlib import Path

import joblib
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_firewall.models.decisions import DecisionThresholds
from genome_firewall.models.metrics import binary_metrics, no_call_metrics
from genome_firewall.models.training import grouped_train_test_split, load_model_table


def build_candidates(random_state: int) -> dict[str, object]:
    return {
        "logistic_regression": Pipeline(
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
        ),
        "linear_svc": Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "clf",
                    LinearSVC(
                        class_weight="balanced",
                        max_iter=5000,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced_subsample",
            max_features="sqrt",
            min_samples_leaf=2,
            random_state=random_state,
            n_jobs=-1,
        ),
        "extra_trees": ExtraTreesClassifier(
            n_estimators=400,
            class_weight="balanced",
            max_features="sqrt",
            min_samples_leaf=2,
            random_state=random_state,
            n_jobs=-1,
        ),
        "gradient_boosting": GradientBoostingClassifier(
            n_estimators=150,
            learning_rate=0.05,
            max_depth=2,
            random_state=random_state,
        ),
    }


def calibrated(estimator: object) -> CalibratedClassifierCV:
    return CalibratedClassifierCV(estimator=estimator, method="sigmoid", cv=3)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare small-data classical models on grouped AMR splits.")
    parser.add_argument("--labels", type=Path, default=ROOT / "data" / "processed" / "ecoli_cohort_labels.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data" / "processed" / "ecoli_genome_metadata.csv")
    parser.add_argument(
        "--features",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refined_features.csv",
    )
    parser.add_argument("--feature-prefix", default="refined_")
    parser.add_argument("--out", type=Path, default=ROOT / "reports" / "metrics" / "classical_model_comparison.csv")
    parser.add_argument("--best-out", type=Path, default=ROOT / "reports" / "metrics" / "best_classical_models.csv")
    parser.add_argument("--model-dir", type=Path, default=ROOT / "models" / "classical_comparison")
    parser.add_argument("--min-class-count", type=int, default=8)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--likely-to-work-max", type=float, default=0.30)
    parser.add_argument("--likely-to-fail-min", type=float, default=0.70)
    args = parser.parse_args()

    thresholds = DecisionThresholds(
        likely_to_work_max=args.likely_to_work_max,
        likely_to_fail_min=args.likely_to_fail_min,
    )
    table, feature_columns = load_model_table(args.labels, args.metadata, args.features, args.feature_prefix)
    candidates = build_candidates(args.random_state)
    rows: list[dict] = []

    args.model_dir.mkdir(parents=True, exist_ok=True)
    for antibiotic, antibiotic_rows in table.groupby("antibiotic"):
        counts = antibiotic_rows["y"].value_counts().to_dict()
        if len(counts) < 2 or min(counts.values()) < args.min_class_count:
            print(f"Skipping {antibiotic}: insufficient class balance {counts}", flush=True)
            continue

        train, test = grouped_train_test_split(antibiotic_rows, args.random_state, args.test_size)
        if train["y"].nunique() < 2 or test["y"].nunique() < 2:
            print(f"Skipping {antibiotic}: grouped split lost one class", flush=True)
            continue

        for model_name, estimator in candidates.items():
            try:
                model = calibrated(estimator)
                model.fit(train[feature_columns], train["y"])
                probabilities = model.predict_proba(test[feature_columns])[:, 1]
            except Exception as exc:
                rows.append(
                    {
                        "antibiotic": antibiotic,
                        "model": model_name,
                        "status": "failed",
                        "error": str(exc),
                    }
                )
                print(f"Failed {antibiotic} / {model_name}: {exc}", flush=True)
                continue

            metrics = binary_metrics(test["y"].to_numpy(), probabilities)
            metrics.update(no_call_metrics(test["y"].to_numpy(), probabilities, thresholds))
            row = {
                "antibiotic": antibiotic,
                "model": model_name,
                "status": "trained",
                "error": "",
                "train_rows": len(train),
                "test_rows": len(test),
                "train_groups": train["group"].nunique(),
                "test_groups": test["group"].nunique(),
                "overlapping_groups": len(set(train["group"].astype(str)) & set(test["group"].astype(str))),
                **metrics,
            }
            rows.append(row)
            model_path = args.model_dir / f"{antibiotic.replace('/', '_')}__{model_name}.joblib"
            joblib.dump(model, model_path)
            print(f"Trained {antibiotic} / {model_name}", flush=True)

    comparison = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.out, index=False)

    trained = comparison[comparison["status"].eq("trained")].copy()
    if trained.empty:
        best = pd.DataFrame()
    else:
        best = (
            trained.sort_values(
                ["antibiotic", "balanced_accuracy", "called_accuracy", "brier_score"],
                ascending=[True, False, False, True],
            )
            .groupby("antibiotic", as_index=False)
            .head(1)
        )
    best.to_csv(args.best_out, index=False)

    print(f"Wrote {args.out}", flush=True)
    print(f"Wrote {args.best_out}", flush=True)
    if not best.empty:
        print(best[["antibiotic", "model", "balanced_accuracy", "brier_score", "no_call_rate"]].to_string(index=False))


if __name__ == "__main__":
    main()

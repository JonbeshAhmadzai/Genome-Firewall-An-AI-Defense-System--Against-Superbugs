from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
LABELS_PATH = ROOT / "data" / "processed" / "ecoli_cohort_labels.csv"
METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
FEATURES_PATH = ROOT / "data" / "interim" / "ecoli_kmer_features.csv"
METRICS_PATH = ROOT / "reports" / "metrics" / "kmer_baseline_metrics.csv"
MODELS_DIR = ROOT / "models" / "kmer_baseline"
MANIFEST_PATH = MODELS_DIR / "manifest.json"


def score_model(y_true: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = (probabilities >= 0.5).astype(int)
    metrics = {
        "balanced_accuracy": balanced_accuracy_score(y_true, predictions),
        "resistant_recall": recall_score(y_true, predictions, pos_label=1, zero_division=0),
        "susceptible_recall": recall_score(y_true, predictions, pos_label=0, zero_division=0),
        "f1": f1_score(y_true, predictions, zero_division=0),
        "brier_score": brier_score_loss(y_true, probabilities),
    }
    if len(set(y_true)) == 2:
        metrics["auroc"] = roc_auc_score(y_true, probabilities)
        metrics["pr_auc"] = average_precision_score(y_true, probabilities)
    else:
        metrics["auroc"] = np.nan
        metrics["pr_auc"] = np.nan
    return metrics


def main() -> None:
    labels = pd.read_csv(LABELS_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    features = pd.read_csv(FEATURES_PATH)

    labels["genome_id"] = labels["genome_id"].astype(str)
    metadata["genome_id"] = metadata["genome_id"].astype(str)
    features["genome_id"] = features["genome_id"].astype(str)

    feature_columns = [column for column in features.columns if column.startswith("kmer_")]
    merged = labels.merge(features, on="genome_id", how="inner").merge(
        metadata[["genome_id", "cgmlst_hc100"]], on="genome_id", how="left"
    ).copy()
    merged["y"] = merged["label"].map({"susceptible": 0, "resistant": 1})
    merged["group"] = merged["cgmlst_hc100"].fillna(merged["genome_id"]).astype(str)

    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)

    metric_rows = []
    manifest = {
        "feature_type": "temporary_4mer_frequency",
        "warning": "Plumbing baseline only. Replace with AMRFinderPlus gene/mutation features for the challenge submission.",
        "models": {},
    }

    for antibiotic, group in merged.groupby("antibiotic"):
        counts = group["y"].value_counts()
        if len(counts) < 2 or counts.min() < 8:
            print(f"Skipping {antibiotic}: insufficient class balance {counts.to_dict()}")
            continue

        splitter = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=42)
        train_idx, test_idx = next(splitter.split(group, group["y"], groups=group["group"]))
        train = group.iloc[train_idx]
        test = group.iloc[test_idx]

        if train["y"].nunique() < 2 or test["y"].nunique() < 2:
            print(f"Skipping {antibiotic}: grouped split lost a class")
            continue

        model = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        class_weight="balanced",
                        max_iter=2000,
                        solver="liblinear",
                    ),
                ),
            ]
        )
        model.fit(train[feature_columns], train["y"])
        probabilities = model.predict_proba(test[feature_columns])[:, 1]
        metrics = score_model(test["y"].to_numpy(), probabilities)
        metrics.update(
            {
                "antibiotic": antibiotic,
                "train_rows": len(train),
                "test_rows": len(test),
                "train_groups": train["group"].nunique(),
                "test_groups": test["group"].nunique(),
            }
        )
        metric_rows.append(metrics)

        model_path = MODELS_DIR / f"{antibiotic.replace('/', '_')}.joblib"
        joblib.dump(model, model_path)
        manifest["models"][antibiotic] = {
            "path": str(model_path.relative_to(ROOT)),
            "feature_columns": feature_columns,
        }
        print(f"Trained {antibiotic}")

    pd.DataFrame(metric_rows).to_csv(METRICS_PATH, index=False)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {METRICS_PATH}")
    print(f"Wrote {MANIFEST_PATH}")


if __name__ == "__main__":
    main()

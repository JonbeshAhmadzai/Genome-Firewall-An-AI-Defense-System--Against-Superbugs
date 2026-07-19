from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from genome_firewall.models.decisions import DecisionThresholds, classify_probability, evidence_category
from genome_firewall.models.metrics import binary_metrics, no_call_metrics


LABEL_MAP = {"susceptible": 0, "resistant": 1}
SAFETY_WARNING = "Research prototype only. Confirm with standard laboratory testing."


def load_model_table(
    labels_path: Path,
    metadata_path: Path,
    features_path: Path,
    feature_prefix: str,
) -> tuple[pd.DataFrame, list[str]]:
    labels = pd.read_csv(labels_path)
    metadata = pd.read_csv(metadata_path)
    features = pd.read_csv(features_path)

    for frame in (labels, metadata, features):
        frame["genome_id"] = frame["genome_id"].astype(str)

    feature_columns = [column for column in features.columns if column.startswith(feature_prefix)]
    if not feature_columns:
        raise ValueError(f"No feature columns found with prefix {feature_prefix!r}")

    group_columns = ["genome_id"]
    if "cgmlst_hc100" in metadata.columns:
        group_columns.append("cgmlst_hc100")
    quality_columns = [
        column
        for column in ["genome_quality", "checkm_completeness", "checkm_contamination"]
        if column in metadata.columns
    ]

    merged = labels.merge(features, on="genome_id", how="inner").merge(
        metadata[group_columns + quality_columns], on="genome_id", how="left"
    ).copy()
    merged["y"] = merged["label"].map(LABEL_MAP)
    merged = merged.dropna(subset=["y"]).copy()
    merged["y"] = merged["y"].astype(int)
    if "cgmlst_hc100" in merged.columns:
        merged["group"] = merged["cgmlst_hc100"].fillna(merged["genome_id"]).astype(str)
    else:
        merged["group"] = merged["genome_id"]
    return merged, feature_columns


def build_classifier(random_state: int) -> CalibratedClassifierCV:
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
    return CalibratedClassifierCV(base, method="sigmoid", cv=3)


def grouped_train_test_split(group: pd.DataFrame, random_state: int, test_size: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(splitter.split(group, group["y"], groups=group["group"]))
    return group.iloc[train_idx].copy(), group.iloc[test_idx].copy()


def train_models(
    table: pd.DataFrame,
    feature_columns: list[str],
    output_dir: Path,
    feature_type: str,
    min_class_count: int = 8,
    random_state: int = 42,
    test_size: float = 0.25,
    thresholds: DecisionThresholds | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = thresholds or DecisionThresholds()
    metric_rows = []
    prediction_rows = []
    split_audit_rows = []
    manifest = {
        "feature_type": feature_type,
        "thresholds": asdict(thresholds),
        "label_map": LABEL_MAP,
        "safety_warning": SAFETY_WARNING,
        "models": {},
        "skipped": {},
    }

    for antibiotic, antibiotic_rows in table.groupby("antibiotic"):
        counts = antibiotic_rows["y"].value_counts().to_dict()
        if len(counts) < 2 or min(counts.values()) < min_class_count:
            manifest["skipped"][antibiotic] = f"insufficient class balance: {counts}"
            print(f"Skipping {antibiotic}: {manifest['skipped'][antibiotic]}")
            continue

        train, test = grouped_train_test_split(antibiotic_rows, random_state, test_size)
        overlap = sorted(set(train["group"].astype(str)) & set(test["group"].astype(str)))
        if train["y"].nunique() < 2 or test["y"].nunique() < 2:
            manifest["skipped"][antibiotic] = "grouped split lost one class"
            print(f"Skipping {antibiotic}: {manifest['skipped'][antibiotic]}")
            continue
        train_counts = train["y"].value_counts().to_dict()
        if min(train_counts.values()) < 3:
            manifest["skipped"][antibiotic] = f"not enough training samples per class for calibration: {train_counts}"
            print(f"Skipping {antibiotic}: {manifest['skipped'][antibiotic]}")
            continue

        model = build_classifier(random_state)
        model.fit(train[feature_columns], train["y"])
        probabilities = model.predict_proba(test[feature_columns])[:, 1]

        metrics = binary_metrics(test["y"].to_numpy(), probabilities)
        metrics.update(no_call_metrics(test["y"].to_numpy(), probabilities, thresholds))
        metrics.update(
            {
                "antibiotic": antibiotic,
                "train_rows": len(train),
                "test_rows": len(test),
                "train_groups": train["group"].nunique(),
                "test_groups": test["group"].nunique(),
                "overlapping_groups": len(overlap),
                "resistant_train": int((train["y"] == 1).sum()),
                "susceptible_train": int((train["y"] == 0).sum()),
                "resistant_test": int((test["y"] == 1).sum()),
                "susceptible_test": int((test["y"] == 0).sum()),
            }
        )
        metric_rows.append(metrics)
        split_audit_rows.append(
            {
                "antibiotic": antibiotic,
                "group_column": "cgmlst_hc100" if "cgmlst_hc100" in antibiotic_rows.columns else "genome_id",
                "train_rows": len(train),
                "test_rows": len(test),
                "train_genomes": train["genome_id"].nunique(),
                "test_genomes": test["genome_id"].nunique(),
                "train_groups": int(train["group"].nunique()),
                "test_groups": int(test["group"].nunique()),
                "overlapping_groups": len(overlap),
                "overlap_examples": ", ".join(overlap[:10]),
            }
        )

        for row, probability in zip(test.itertuples(index=False), probabilities):
            decision, confidence = classify_probability(float(probability), thresholds)
            prediction_rows.append(
                {
                    "genome_id": row.genome_id,
                    "antibiotic": antibiotic,
                    "true_label": row.label,
                    "resistant_probability": float(probability),
                    "decision": decision,
                    "confidence": confidence,
                    "evidence_category": evidence_category(feature_type),
                    "safety_warning": SAFETY_WARNING,
                }
            )

        model_path = output_dir / f"{antibiotic.replace('/', '_')}.joblib"
        joblib.dump(model, model_path)
        manifest["models"][antibiotic] = {
            "path": str(model_path),
            "feature_columns": feature_columns,
            "train_rows": len(train),
            "test_rows": len(test),
            "train_groups": int(train["group"].nunique()),
            "test_groups": int(test["group"].nunique()),
        }
        print(f"Trained {antibiotic}")

    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), manifest, pd.DataFrame(split_audit_rows)


def write_manifest(manifest: dict, path: Path, root: Path) -> None:
    serializable = json.loads(json.dumps(manifest))
    for model in serializable["models"].values():
        model["path"] = str(Path(model["path"]).resolve().relative_to(root.resolve()))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

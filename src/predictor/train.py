"""Train one calibrated logistic model per validated species/antibiotic pair."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import PredefinedSplit

from .grouped_split import SplitAssignment, grouped_three_way_split, split_receipt


def _calibrator(base: LogisticRegression, calibration_rows: int) -> CalibratedClassifierCV:
    """Create a calibration wrapper compatible with supported sklearn versions."""

    try:
        from sklearn.frozen import FrozenEstimator

        # One predefined fold makes every calibration row a held-out
        # calibration observation for the already-fitted base estimator. This
        # avoids an implicit 5-fold split, which is invalid for the MVP-sized
        # calibration sets.
        cv = PredefinedSplit(np.zeros(calibration_rows, dtype=int))
        return CalibratedClassifierCV(FrozenEstimator(base), method="sigmoid", cv=cv)
    except ImportError:  # scikit-learn 1.4/1.5 compatibility
        return CalibratedClassifierCV(base, method="sigmoid", cv="prefit")


def prepare_target_frame(
    labels: pd.DataFrame,
    features: pd.DataFrame,
    groups: pd.DataFrame,
    *,
    species: str,
    antibiotic: str,
) -> pd.DataFrame:
    """Join labels, features, and homology groups for one target pair."""

    required_labels = {"genome_id", "species", "antibiotic", "phenotype"}
    missing = required_labels - set(labels.columns)
    if missing:
        raise ValueError(f"Labels are missing columns: {sorted(missing)}")
    if "genome_id" not in features.columns or "genome_id" not in groups.columns:
        raise ValueError("Features and groups must both contain genome_id")
    selected = labels[labels["species"].eq(species) & labels["antibiotic"].eq(antibiotic)].copy()
    selected["y"] = selected["phenotype"].map({"susceptible": 0, "resistant": 1})
    selected = selected[selected["y"].notna()].copy()
    merged = selected.merge(features, on="genome_id", how="inner", validate="one_to_one")
    merged = merged.merge(groups[["genome_id", "homology_group_id"]], on="genome_id", how="inner", validate="one_to_one")
    if merged.empty:
        raise ValueError(f"No aligned rows for {species} / {antibiotic}")
    return merged


def train_target(
    target_frame: pd.DataFrame,
    *,
    output_dir: Path,
    species: str,
    antibiotic: str,
    random_state: int = 42,
) -> dict[str, object]:
    """Fit, calibrate, and save a single target model with a split receipt."""

    assignment: SplitAssignment = grouped_three_way_split(target_frame, random_state=random_state)
    split = assignment.frame
    excluded = {"genome_id", "species", "antibiotic", "phenotype", "y", "homology_group_id", "split", "_original_index"}
    feature_columns = [
        column
        for column in split.columns
        if column not in excluded and pd.api.types.is_numeric_dtype(split[column])
    ]
    if not feature_columns:
        raise ValueError("No numeric feature columns available for model training")
    X = split[feature_columns].apply(pd.to_numeric, errors="raise")
    y = split["y"].astype(int)
    train_mask = split["split"].eq("train")
    calibration_mask = split["split"].eq("calibration")
    if y[train_mask].nunique() < 2 or y[calibration_mask].nunique() < 2:
        raise ValueError("Training and calibration groups must both contain Resistant and Susceptible classes")

    base = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=random_state)
    base.fit(X.loc[train_mask], y.loc[train_mask])
    calibrated = _calibrator(base, int(calibration_mask.sum()))
    calibrated.fit(X.loc[calibration_mask], y.loc[calibration_mask])

    output_dir.mkdir(parents=True, exist_ok=True)
    slug = f"{species.lower().replace(' ', '_')}__{antibiotic.lower().replace(' ', '_')}"
    model_path = output_dir / f"{slug}.joblib"
    receipt_path = output_dir / f"{slug}.json"
    split_path = output_dir / f"{slug}__splits.csv"
    joblib.dump({"model": calibrated, "feature_columns": feature_columns, "species": species, "antibiotic": antibiotic}, model_path)
    split.to_csv(split_path, index=False)
    receipt = {
        "species": species,
        "antibiotic": antibiotic,
        "feature_columns": feature_columns,
        "split": split_receipt(assignment),
        "model_path": str(model_path),
        "split_path": str(split_path),
        "classes": {"susceptible": 0, "resistant": 1},
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt

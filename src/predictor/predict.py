"""Inference helpers for an approved calibrated target model."""

from __future__ import annotations

from pathlib import Path
from collections.abc import Iterable

import joblib
import numpy as np
import pandas as pd

from src.trust.calibration import apply_no_call
from .target_gate import TargetGateDecision, apply_target_gate


def load_model(path: Path) -> dict[str, object]:
    artifact = joblib.load(path)
    if not isinstance(artifact, dict) or "model" not in artifact or "feature_columns" not in artifact:
        raise ValueError(f"Invalid model artifact: {path}")
    return artifact


def predict_target(
    model_artifact: dict[str, object],
    features: pd.DataFrame,
    *,
    target_present: bool | None | Iterable[bool | None],
    target_names: tuple[str, ...] = (),
    min_confidence: float = 0.70,
) -> pd.DataFrame:
    """Predict from one calibrated model and apply trust + target gates."""

    feature_columns = list(model_artifact["feature_columns"])
    missing = [column for column in feature_columns if column not in features.columns]
    if missing:
        raise ValueError(f"Input feature matrix is missing model columns: {missing[:5]}")
    model = model_artifact["model"]
    probabilities = np.asarray(model.predict_proba(features[feature_columns]), dtype=float)[:, 1]
    decisions = apply_no_call(probabilities, min_confidence=min_confidence)
    if isinstance(target_present, (bool, type(None))):
        target_values = [target_present] * len(decisions)
    else:
        target_values = list(target_present)
        if len(target_values) != len(decisions):
            raise ValueError("target_present iterable must match the number of feature rows")
    gated: list[TargetGateDecision] = [
        apply_target_gate(prediction, value, target_names=target_names)
        for prediction, value in zip(decisions["prediction"], target_values)
    ]
    result = decisions.copy()
    result["prediction_before_target_gate"] = result["prediction"]
    result["prediction"] = [decision.prediction for decision in gated]
    result["target_status"] = [decision.target_status for decision in gated]
    result["target_gate_reason"] = [decision.reason for decision in gated]
    if "genome_id" in features.columns:
        result.insert(0, "genome_id", features["genome_id"].to_numpy())
    return result

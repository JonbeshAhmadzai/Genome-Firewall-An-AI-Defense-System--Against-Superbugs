"""Challenge-aligned held-out metrics for one antibiotic model."""

from __future__ import annotations

from typing import Iterable

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    recall_score,
    roc_auc_score,
)

from src.trust.calibration import apply_no_call, reliability_table


def scorecard(
    y_true: Iterable[int],
    probability_resistant: Iterable[float],
    *,
    min_confidence: float = 0.70,
) -> dict[str, float | int]:
    truth = np.asarray(list(y_true), dtype=int)
    probabilities = np.asarray(list(probability_resistant), dtype=float)
    if len(truth) == 0 or len(truth) != len(probabilities):
        raise ValueError("Truth and probability arrays must be non-empty and equal length")
    prediction = (probabilities >= 0.5).astype(int)
    decisions = apply_no_call(probabilities, min_confidence=min_confidence)
    called = ~decisions["no_call"].to_numpy()
    reliability, brier = reliability_table(truth, probabilities)
    metrics: dict[str, float | int] = {
        "n_test": int(len(truth)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, prediction)),
        "resistant_recall": float(recall_score(truth, prediction, pos_label=1, zero_division=0)),
        "susceptible_recall": float(recall_score(truth, prediction, pos_label=0, zero_division=0)),
        "f1": float(f1_score(truth, prediction, zero_division=0)),
        "auroc": float(roc_auc_score(truth, probabilities)) if np.unique(truth).size == 2 else float("nan"),
        "pr_auc": float(average_precision_score(truth, probabilities)) if np.unique(truth).size == 2 else float("nan"),
        "brier": brier,
        "no_call_rate": float((~called).mean()),
        "called_accuracy": float((prediction[called] == truth[called]).mean()) if called.any() else float("nan"),
        "reliability_bins": int(len(reliability)),
    }
    return metrics

from __future__ import annotations

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    recall_score,
    roc_auc_score,
)


def binary_metrics(y_true: np.ndarray, resistant_probability: np.ndarray) -> dict[str, float]:
    predictions = (resistant_probability >= 0.5).astype(int)
    metrics = {
        "balanced_accuracy": balanced_accuracy_score(y_true, predictions),
        "resistant_recall": recall_score(y_true, predictions, pos_label=1, zero_division=0),
        "susceptible_recall": recall_score(y_true, predictions, pos_label=0, zero_division=0),
        "f1": f1_score(y_true, predictions, zero_division=0),
        "brier_score": brier_score_loss(y_true, resistant_probability),
    }
    if len(set(y_true)) == 2:
        metrics["auroc"] = roc_auc_score(y_true, resistant_probability)
        metrics["pr_auc"] = average_precision_score(y_true, resistant_probability)
    else:
        metrics["auroc"] = np.nan
        metrics["pr_auc"] = np.nan
    return metrics


def no_call_metrics(y_true: np.ndarray, resistant_probability: np.ndarray, thresholds=None) -> dict[str, float]:
    from genome_firewall.models.decisions import classify_probability

    decisions = np.array([classify_probability(probability, thresholds)[0] for probability in resistant_probability])
    called = decisions != "no_call"
    result = {
        "no_call_rate": float(1.0 - called.mean()) if len(called) else np.nan,
        "called_rows": int(called.sum()),
    }
    if called.any():
        called_predictions = np.where(decisions[called] == "likely_to_fail", 1, 0)
        result["called_accuracy"] = float((called_predictions == y_true[called]).mean())
    else:
        result["called_accuracy"] = np.nan
    return result

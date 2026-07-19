"""Calibration, abstention, and reliability calculations."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss


def confidence_from_probability(probability_resistant: Iterable[float]) -> np.ndarray:
    probabilities = np.asarray(list(probability_resistant), dtype=float)
    if np.any((probabilities < 0) | (probabilities > 1)):
        raise ValueError("Probabilities must be between 0 and 1")
    return np.maximum(probabilities, 1.0 - probabilities)


def apply_no_call(
    probability_resistant: Iterable[float],
    *,
    min_confidence: float = 0.70,
    out_of_distribution: Iterable[bool] | None = None,
) -> pd.DataFrame:
    """Return a prediction table with conservative no-call decisions."""

    probabilities = np.asarray(list(probability_resistant), dtype=float)
    confidence = confidence_from_probability(probabilities)
    if not 0.5 <= min_confidence <= 1:
        raise ValueError("min_confidence must be between 0.5 and 1")
    if out_of_distribution is None:
        ood = np.zeros(len(probabilities), dtype=bool)
    else:
        ood = np.asarray(list(out_of_distribution), dtype=bool)
        if len(ood) != len(probabilities):
            raise ValueError("out_of_distribution must match the probability length")
    called = (confidence >= min_confidence) & ~ood
    labels = np.where(probabilities >= 0.5, "likely to fail", "likely to work")
    labels = np.where(called, labels, "no-call")
    return pd.DataFrame(
        {
            "probability_resistant": probabilities,
            "confidence": confidence,
            "prediction": labels,
            "no_call": ~called,
            "out_of_distribution": ood,
        }
    )


def reliability_table(
    y_true: Iterable[int],
    probability_resistant: Iterable[float],
    *,
    bins: int = 10,
) -> tuple[pd.DataFrame, float]:
    """Return equal-width reliability bins and the Brier score."""

    truth = np.asarray(list(y_true), dtype=int)
    probabilities = np.asarray(list(probability_resistant), dtype=float)
    if len(truth) != len(probabilities) or len(truth) == 0:
        raise ValueError("Truth and probability arrays must be non-empty and equal length")
    if bins < 2:
        raise ValueError("At least two reliability bins are required")
    edges = np.linspace(0, 1, bins + 1)
    bucket = np.clip(np.digitize(probabilities, edges[1:-1], right=False), 0, bins - 1)
    rows = []
    for index in range(bins):
        mask = bucket == index
        rows.append(
            {
                "bin": index,
                "lower": edges[index],
                "upper": edges[index + 1],
                "count": int(mask.sum()),
                "mean_confidence": float(probabilities[mask].mean()) if mask.any() else np.nan,
                "observed_resistant_rate": float(truth[mask].mean()) if mask.any() else np.nan,
            }
        )
    return pd.DataFrame(rows), float(brier_score_loss(truth, probabilities))

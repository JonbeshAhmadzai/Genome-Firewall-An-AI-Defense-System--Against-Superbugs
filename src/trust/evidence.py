"""Evidence typing that keeps curated biology separate from model association."""

from __future__ import annotations

from typing import Iterable


def classify_evidence(
    observed_features: Iterable[str],
    known_features: Iterable[str],
) -> str:
    """Classify evidence as known determinant, statistical association, or none."""

    observed = {str(value).strip() for value in observed_features if str(value).strip()}
    known = {str(value).strip() for value in known_features if str(value).strip()}
    if observed & known:
        return "known_resistance_determinant"
    if observed:
        return "statistical_association"
    return "no_known_signal"

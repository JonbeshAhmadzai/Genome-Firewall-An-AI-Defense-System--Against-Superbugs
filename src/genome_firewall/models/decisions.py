from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DecisionThresholds:
    likely_to_work_max: float = 0.30
    likely_to_fail_min: float = 0.70


def classify_probability(
    resistant_probability: float,
    thresholds: DecisionThresholds = DecisionThresholds(),
) -> tuple[str, float]:
    """Map resistant probability to a challenge-facing decision and confidence."""
    if resistant_probability >= thresholds.likely_to_fail_min:
        return "likely_to_fail", resistant_probability
    if resistant_probability <= thresholds.likely_to_work_max:
        return "likely_to_work", 1.0 - resistant_probability
    return "no_call", max(resistant_probability, 1.0 - resistant_probability)


def evidence_category(feature_type: str) -> str:
    if "amrfinder" in feature_type.lower():
        return "known resistance gene or DNA change detected when supporting features are present"
    return "statistical association only"


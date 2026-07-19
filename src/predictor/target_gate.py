"""Conservative molecular-target compatibility gate.

This module does not predict phenotype. It only prevents a model from saying
"likely to work" when the drug's molecular target is absent or unverified.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class TargetGateDecision:
    prediction: str
    target_status: str
    reason: str


def apply_target_gate(
    prediction: str,
    target_present: bool | None,
    *,
    target_names: tuple[str, ...] = (),
) -> TargetGateDecision:
    """Apply the safety rule to a model prediction.

    ``target_present`` is intentionally tri-state. ``None`` means that the
    feature/sequence check was inconclusive, which is treated like absence for
    a proposed susceptible/work result.
    """

    normalized = prediction.strip().lower()
    allowed = {"likely to work", "likely to fail", "no-call"}
    if normalized not in allowed:
        raise ValueError(f"Unsupported prediction {prediction!r}; expected one of {sorted(allowed)}")

    target_text = ", ".join(target_names) if target_names else "the configured molecular target"
    if normalized == "likely to work" and target_present is not True:
        status = "absent" if target_present is False else "uncertain"
        return TargetGateDecision(
            prediction="no-call",
            target_status=status,
            reason=f"Cannot call likely to work: {target_text} was not verified as present.",
        )

    status = "present" if target_present is True else "absent" if target_present is False else "uncertain"
    if normalized == "no-call":
        return TargetGateDecision(prediction="no-call", target_status=status, reason="Model evidence was not strong enough.")
    if target_present is False:
        return TargetGateDecision(prediction="no-call", target_status=status, reason=f"Target compatibility failed for {target_text}.")
    if target_present is None:
        return TargetGateDecision(prediction="no-call", target_status=status, reason=f"Target compatibility is uncertain for {target_text}.")
    return TargetGateDecision(prediction=normalized, target_status=status, reason="Target compatibility verified.")


def target_presence_from_features(
    feature_names: Mapping[str, bool],
    target_names: tuple[str, ...],
) -> bool | None:
    """Resolve target presence from normalized annotation names.

    A positive match is sufficient for this conservative feature-level check.
    An empty feature map is unknown rather than absent, because an annotation
    failure must not be mistaken for biological absence.
    """

    if not feature_names:
        return None
    normalized = {name.strip().lower(): bool(value) for name, value in feature_names.items()}
    for target in target_names:
        if normalized.get(target.strip().lower()) is True:
            return True
    return False

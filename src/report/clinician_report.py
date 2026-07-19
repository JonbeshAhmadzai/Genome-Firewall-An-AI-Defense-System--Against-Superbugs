"""Structured, non-decisional report text for the eventual application layer."""

from __future__ import annotations

import pandas as pd


LAB_DISCLAIMER = "This is decision support — confirm with standard lab testing."


def build_report(predictions: pd.DataFrame, *, species: str, antibiotic: str) -> str:
    """Render model output without adding a label or biological claim."""

    required = {"prediction", "confidence", "target_status"}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"Predictions are missing report fields: {sorted(missing)}")
    lines = [
        f"Genome Firewall report for {species} — {antibiotic}",
        "",
        "This report summarizes a calibrated research model. It does not make a treatment decision.",
    ]
    for _, row in predictions.iterrows():
        genome = row.get("genome_id", "sample")
        lines.append(
            f"{genome}: {row['prediction']} (confidence {float(row['confidence']):.2f}; "
            f"molecular target status: {row['target_status']})."
        )
    lines.extend(["", LAB_DISCLAIMER])
    return "\n".join(lines)

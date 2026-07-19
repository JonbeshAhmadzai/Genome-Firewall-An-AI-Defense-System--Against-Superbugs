"""Evaluate a saved target model on its untouched grouped test split."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.trust.calibration import apply_no_call, reliability_table
from .scorecard import scorecard


def evaluate_artifact(
    model_artifact: dict[str, object],
    split_frame: pd.DataFrame,
    *,
    output_dir: Path | None = None,
    min_confidence: float = 0.70,
) -> dict[str, object]:
    """Score one model and return metrics plus per-genetic-group results."""

    if "split" not in split_frame.columns:
        raise ValueError("Evaluation frame must contain split assignments")
    test = split_frame[split_frame["split"].eq("test")].copy()
    feature_columns = list(model_artifact["feature_columns"])
    missing = [column for column in feature_columns if column not in test.columns]
    if missing:
        raise ValueError(f"Evaluation frame is missing model features: {missing[:5]}")
    if "y" not in test.columns:
        raise ValueError("Evaluation frame must contain numeric y labels")
    probabilities = model_artifact["model"].predict_proba(test[feature_columns])[:, 1]
    metrics = scorecard(test["y"].astype(int), probabilities, min_confidence=min_confidence)
    reliability, brier = reliability_table(test["y"].astype(int), probabilities)
    decisions = apply_no_call(probabilities, min_confidence=min_confidence)
    test["probability_resistant"] = probabilities
    test["prediction"] = decisions["prediction"].to_numpy()
    test["confidence"] = decisions["confidence"].to_numpy()
    groups = []
    if "homology_group_id" in test.columns:
        for group, subset in test.groupby("homology_group_id", sort=True):
            called = subset["prediction"].ne("no-call")
            observed = subset["phenotype"].map({"susceptible": 0, "resistant": 1}) if "phenotype" in subset.columns else None
            predicted = subset["prediction"].eq("likely to fail").astype(int)
            groups.append(
                {
                    "homology_group_id": str(group),
                    "n": int(len(subset)),
                    "called": int(called.sum()),
                    "no_call_rate": float((~called).mean()),
                    "called_accuracy": float((predicted.loc[called].eq(observed.loc[called]).mean()))
                    if called.any() and observed is not None
                    else float("nan"),
                }
            )
    result: dict[str, object] = {
        "metrics": metrics,
        "brier": brier,
        "group_breakdown": groups,
        "test_predictions": test,
    }
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        reliability.to_csv(output_dir / "reliability.csv", index=False)
        test.to_csv(output_dir / "test_predictions.csv", index=False)
        receipt = {"metrics": metrics, "brier": brier, "group_breakdown": groups}
        (output_dir / "validation_report.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n", encoding="utf-8")
    return result

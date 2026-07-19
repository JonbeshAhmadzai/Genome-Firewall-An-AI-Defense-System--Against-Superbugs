from __future__ import annotations

import json
from pathlib import Path

import joblib
import pandas as pd

from genome_firewall.models.decisions import DecisionThresholds, classify_probability, evidence_category
from genome_firewall.models.training import SAFETY_WARNING


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def predict_feature_table(features: pd.DataFrame, manifest_path: Path, root: Path) -> pd.DataFrame:
    manifest = load_manifest(manifest_path)
    thresholds = DecisionThresholds(**manifest.get("thresholds", {}))
    feature_type = manifest.get("feature_type", "unknown")
    rows = []
    features = features.copy()
    features["genome_id"] = features["genome_id"].astype(str)

    for antibiotic, model_info in manifest.get("models", {}).items():
        model_path = root / model_info["path"]
        model = joblib.load(model_path)
        feature_columns = model_info["feature_columns"]
        missing = [column for column in feature_columns if column not in features.columns]
        if missing:
            raise ValueError(f"Missing {len(missing)} feature columns for {antibiotic}")

        probabilities = model.predict_proba(features[feature_columns])[:, 1]
        for genome_id, probability in zip(features["genome_id"], probabilities):
            decision, confidence = classify_probability(float(probability), thresholds)
            rows.append(
                {
                    "genome_id": genome_id,
                    "antibiotic": antibiotic,
                    "resistant_probability": float(probability),
                    "decision": decision,
                    "confidence": confidence,
                    "evidence_category": evidence_category(feature_type),
                    "safety_warning": manifest.get("safety_warning", SAFETY_WARNING),
                }
            )
    return pd.DataFrame(rows)


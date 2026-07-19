from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class TargetGateResult:
    compatible: bool
    status: str
    reason: str
    molecular_targets: str
    target_present: bool | None = None
    matched_target_markers: str = ""


def evaluate_target_gate(
    antibiotic: str,
    species: str,
    antibiotic_config: dict,
    target_presence_row: pd.Series | None = None,
) -> TargetGateResult:
    supported_species = str(antibiotic_config.get("supported_species", "")).strip().lower()
    observed_species = str(species or "").strip().lower()
    drug_config = antibiotic_config.get("antibiotics", {}).get(antibiotic, {})
    targets = ", ".join(drug_config.get("primary_targets", []))
    rule = drug_config.get("target_presence_rule", "not_configured")

    if not drug_config:
        return TargetGateResult(
            compatible=False,
            status="no_call",
            reason="antibiotic is not configured for target-gate evaluation",
            molecular_targets="",
        )
    if observed_species != supported_species:
        return TargetGateResult(
            compatible=False,
            status="no_call",
            reason=f"species {species!r} is outside configured scope {antibiotic_config.get('supported_species')!r}",
            molecular_targets=targets,
        )

    if target_presence_row is not None:
        target_status = str(target_presence_row.get("target_status", "")).strip()
        target_present_raw = target_presence_row.get("target_present")
        target_present = str(target_present_raw).strip().lower() in {"true", "1", "yes"}
        matched_markers = str(target_presence_row.get("matched_target_markers", "") or "")
        target_reason = str(target_presence_row.get("target_gate_reason", "") or "")
        row_targets = str(target_presence_row.get("molecular_targets", "") or "")
        if row_targets:
            targets = row_targets

        if target_status == "present" or target_present:
            return TargetGateResult(
                compatible=True,
                status="compatible",
                reason=target_reason or "target marker present in genome annotation",
                molecular_targets=targets,
                target_present=True,
                matched_target_markers=matched_markers,
            )
        if target_status == "absent":
            return TargetGateResult(
                compatible=False,
                status="target_absent",
                reason=target_reason or "configured target marker was not found in genome annotation",
                molecular_targets=targets,
                target_present=False,
                matched_target_markers=matched_markers,
            )
        return TargetGateResult(
            compatible=False,
            status="target_unknown",
            reason=target_reason or "target presence could not be verified from genome annotation",
            molecular_targets=targets,
            target_present=None,
            matched_target_markers=matched_markers,
        )

    if rule == "annotation_markers":
        return TargetGateResult(
            compatible=False,
            status="target_unknown",
            reason="target_presence.csv was not provided; run scripts/fetch_bvbrc_genome_features.py and scripts/build_target_presence.py",
            molecular_targets=targets,
            target_present=None,
        )

    if rule == "assumed_present_for_supported_species":
        return TargetGateResult(
            compatible=True,
            status="compatible",
            reason="target assumed present for the supported species; final version should verify from genome annotation",
            molecular_targets=targets,
            target_present=True,
        )
    return TargetGateResult(
        compatible=False,
        status="no_call",
        reason=f"unsupported target-gate rule: {rule}",
        molecular_targets=targets,
    )


def apply_target_gate(
    predictions: pd.DataFrame,
    metadata: pd.DataFrame,
    antibiotic_config: dict,
    target_presence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    frame = predictions.copy()
    metadata = metadata.copy()
    frame["genome_id"] = frame["genome_id"].astype(str)
    metadata["genome_id"] = metadata["genome_id"].astype(str)

    if "species" not in metadata.columns:
        metadata["species"] = ""
    frame = frame.merge(metadata[["genome_id", "species"]], on="genome_id", how="left")

    target_lookup: dict[tuple[str, str], pd.Series] = {}
    if target_presence is not None and not target_presence.empty:
        target_presence = target_presence.copy()
        target_presence["genome_id"] = target_presence["genome_id"].astype(str)
        target_presence["antibiotic"] = target_presence["antibiotic"].astype(str)
        for _, row in target_presence.iterrows():
            target_lookup[(row["genome_id"], row["antibiotic"])] = row

    statuses = []
    reasons = []
    targets = []
    target_present_values = []
    matched_markers = []
    for row in frame.itertuples(index=False):
        target_row = target_lookup.get((str(row.genome_id), str(row.antibiotic)))
        result = evaluate_target_gate(row.antibiotic, row.species, antibiotic_config, target_row)
        statuses.append(result.status)
        reasons.append(result.reason)
        targets.append(result.molecular_targets)
        target_present_values.append(result.target_present)
        matched_markers.append(result.matched_target_markers)
        mask = frame["genome_id"].eq(row.genome_id) & frame["antibiotic"].eq(row.antibiotic)
        if result.status == "no_call":
            frame.loc[mask, "decision"] = "no_call"
            frame.loc[mask, "confidence"] = 0.0
        elif not result.compatible and row.decision == "likely_to_work":
            frame.loc[mask, "decision"] = "no_call"
            frame.loc[mask, "confidence"] = 0.0

    frame["target_gate_status"] = statuses
    frame["target_gate_reason"] = reasons
    frame["molecular_targets"] = targets
    frame["target_present"] = target_present_values
    frame["matched_target_markers"] = matched_markers
    return frame

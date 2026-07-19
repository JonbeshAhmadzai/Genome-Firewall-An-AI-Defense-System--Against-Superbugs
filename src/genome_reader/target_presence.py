"""Verify drug-target compatibility from general genome annotations."""

from __future__ import annotations

import re
from collections.abc import Mapping

import pandas as pd


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip().lower())


def _gene(value: object) -> str:
    return re.sub(r"[^a-z0-9_/-]+", "", _text(value))


def find_target_markers(
    annotations: pd.DataFrame,
    marker_config: Mapping[str, tuple[str, ...]],
) -> list[str]:
    """Return deterministic gene/product matches for one target definition."""

    if annotations.empty:
        return []
    genes = {_gene(value) for value in marker_config.get("genes", ()) if _gene(value)}
    keywords = [_text(value) for value in marker_config.get("product_keywords", ()) if _text(value)]
    observed_genes = {
        _gene(value)
        for value in annotations.get("gene", pd.Series(dtype=str)).dropna()
        if _gene(value)
    }
    products = annotations.get("product", pd.Series(dtype=str)).map(_text)
    matches = {f"gene:{gene}" for gene in genes & observed_genes}
    for keyword in keywords:
        if products.str.contains(re.escape(keyword), regex=True, na=False).any():
            matches.add(f"product:{keyword}")
    return sorted(matches)


def build_target_presence(
    genome_ids: list[str],
    annotations_by_genome: Mapping[str, pd.DataFrame],
    marker_configs: Mapping[str, Mapping[str, tuple[str, ...]]],
    target_names: Mapping[str, tuple[str, ...]],
) -> pd.DataFrame:
    """Build a tri-state genome × antibiotic target-compatibility table."""

    rows: list[dict[str, object]] = []
    for genome_id in genome_ids:
        annotations = annotations_by_genome.get(genome_id)
        for antibiotic, marker_config in marker_configs.items():
            if annotations is None:
                status = "unknown"
                matches: list[str] = []
                reason = "general genome annotation was not available"
                present: bool | None = None
            else:
                matches = find_target_markers(annotations, marker_config)
                present = bool(matches)
                status = "present" if present else "absent"
                reason = (
                    "configured target marker found in general genome annotation"
                    if present
                    else "no configured target marker found in available general genome annotation"
                )
            rows.append(
                {
                    "genome_id": genome_id,
                    "antibiotic": antibiotic,
                    "target_present": present,
                    "target_status": status,
                    "matched_target_markers": "; ".join(matches[:12]),
                    "molecular_targets": "; ".join(target_names.get(antibiotic, ())),
                    "target_gate_reason": reason,
                }
            )
    return pd.DataFrame(rows)

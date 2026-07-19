"""Convert AMRFinderPlus tabular outputs into a sparse presence matrix."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd


def _column(frame: pd.DataFrame, *names: str) -> str | None:
    normalized = {re.sub(r"\s+", " ", column.strip().lower()): column for column in frame.columns}
    for name in names:
        key = re.sub(r"\s+", " ", name.strip().lower())
        if key in normalized:
            return normalized[key]
    return None


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def feature_name(row: pd.Series) -> str:
    """Choose a stable feature identifier from an AMRFinder result row."""

    mutation = _text(row.get("Mutation")) or _text(row.get("Point mutation"))
    gene = _text(row.get("Gene symbol")) or _text(row.get("Element symbol"))
    sequence = (
        _text(row.get("Sequence name"))
        or _text(row.get("Element name"))
        or _text(row.get("Closest reference name"))
    )
    element = _text(row.get("Element type")) or _text(row.get("Type"))
    if mutation:
        return f"mutation:{mutation}"
    if gene:
        return f"gene:{gene}"
    if sequence:
        return f"sequence:{sequence}"
    if element:
        return f"element:{element}"
    return "feature:unidentified"


def parse_amrfinder(path: Path, genome_id: str | None = None) -> tuple[str, set[str], pd.DataFrame]:
    """Parse one AMRFinder TSV and return its ID, feature set, and metadata."""

    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    resolved_id = genome_id or path.stem
    evidence_columns = [
        "genome_id",
        "feature",
        "gene_symbol",
        "mutation",
        "sequence_name",
        "element_type",
        "subtype",
        "amr_class",
        "amr_subclass",
        "method",
        "coverage_reference",
        "identity_reference",
    ]
    if frame.empty:
        return resolved_id, set(), pd.DataFrame(columns=evidence_columns)
    type_column = _column(frame, "Type")
    if type_column:
        # ``--plus`` can also emit stress and virulence hits. They remain in
        # the source TSV, but only resistance determinants are model inputs.
        frame = frame[frame[type_column].str.strip().str.upper().eq("AMR")].copy()
    if frame.empty:
        return resolved_id, set(), pd.DataFrame(columns=evidence_columns)
    frame["feature"] = frame.apply(feature_name, axis=1)
    gene_column = _column(frame, "Gene symbol", "Element symbol")
    mutation_column = _column(frame, "Mutation", "Point mutation")
    sequence_column = _column(frame, "Sequence name", "Element name", "Closest reference name")
    subtype_column = _column(frame, "Subtype", "Element subtype")
    class_column = _column(frame, "Class")
    subclass_column = _column(frame, "Subclass")
    method_column = _column(frame, "Method")
    coverage_column = _column(frame, "% Coverage of reference")
    identity_column = _column(frame, "% Identity to reference")
    metadata = pd.DataFrame(
        {
            "genome_id": resolved_id,
            "feature": frame["feature"],
            "gene_symbol": frame[gene_column] if gene_column else "",
            "mutation": frame[mutation_column] if mutation_column else "",
            "sequence_name": frame[sequence_column] if sequence_column else "",
            "element_type": frame[type_column] if type_column else "AMR",
            "subtype": frame[subtype_column] if subtype_column else "",
            "amr_class": frame[class_column] if class_column else "",
            "amr_subclass": frame[subclass_column] if subclass_column else "",
            "method": frame[method_column] if method_column else "",
            "coverage_reference": frame[coverage_column] if coverage_column else "",
            "identity_reference": frame[identity_column] if identity_column else "",
        }
    )[evidence_columns].drop_duplicates().reset_index(drop=True)
    return resolved_id, set(metadata["feature"]), metadata


def build_feature_matrix(paths: Iterable[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a binary genome × feature matrix and evidence metadata."""

    parsed: list[tuple[str, set[str]]] = []
    metadata: list[pd.DataFrame] = []
    for path in paths:
        genome_id, features, details = parse_amrfinder(Path(path))
        parsed.append((genome_id, features))
        metadata.append(details)
    if not parsed:
        return pd.DataFrame(columns=["genome_id"]), pd.DataFrame()
    vocabulary = sorted({feature for _, features in parsed for feature in features})
    matrix = pd.DataFrame(
        [{"genome_id": genome_id, **{feature: int(feature in features) for feature in vocabulary}} for genome_id, features in parsed]
    )
    evidence = pd.concat(metadata, ignore_index=True) if metadata else pd.DataFrame()
    matrix = add_evidence_features(matrix, evidence)
    return matrix, evidence


def add_evidence_features(matrix: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Add compact per-genome evidence counts without replacing determinants.

    The raw binary gene/mutation columns remain the primary interpretable
    features. Counts provide tree models and linear baselines with robust
    summary signals when individual determinants are sparse.
    """

    frame = matrix.copy()
    if evidence.empty or "genome_id" not in evidence.columns:
        return frame
    evidence = evidence.copy()
    evidence["genome_id"] = evidence["genome_id"].astype(str)
    grouped = evidence.groupby("genome_id", dropna=False)
    summary = grouped.agg(
        amr_hit_count=("feature", "count"),
        amr_class_count=("amr_class", "nunique"),
        amr_subclass_count=("amr_subclass", "nunique"),
        amr_gene_count=("gene_symbol", lambda values: int(values.astype(str).str.strip().ne("").sum())),
        amr_point_mutation_count=("subtype", lambda values: int(values.astype(str).str.upper().str.contains("POINT").sum())),
    ).reset_index()
    return frame.merge(summary, on="genome_id", how="left").fillna(
        {column: 0 for column in summary.columns if column != "genome_id"}
    )


def write_feature_artifacts(paths: Iterable[Path], matrix_path: Path, metadata_path: Path) -> None:
    matrix, evidence = build_feature_matrix(paths)
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(matrix_path, index=False)
    evidence.to_csv(metadata_path, index=False)

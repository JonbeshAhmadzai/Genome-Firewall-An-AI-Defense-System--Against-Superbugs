#!/usr/bin/env python3
"""Fetch and audit BV-BRC laboratory AMR phenotype observations.

This command intentionally stops before downloading genomes. It produces the
metadata receipts required for Human Gate #1 so the team can choose viable
species/antibiotic pairs before spending disk space on FASTA files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping
from urllib.parse import quote

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import targets_config as cfg  # noqa: E402


AMR_TABLE_URL = "ftps://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genome_AMR.txt"
AMR_API_URL = "https://www.bv-brc.org/api/genome_amr/"
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "raw" / "bvbrc"
SCHEMA_VERSION = "1.0"

COLUMN_CANDIDATES: Mapping[str, tuple[str, ...]] = {
    "genome_id": ("genome id", "genome_id"),
    "genome_name": ("genome name", "genome_name"),
    "species": ("species", "organism", "taxon name", "taxon_name"),
    "antibiotic": ("antibiotic", "antimicrobial", "drug"),
    "phenotype": ("resistant phenotype", "resistant_phenotype", "phenotype"),
    "evidence": ("evidence", "amr evidence", "amr_evidence"),
    "measurement": ("measurement",),
    "measurement_sign": ("measurement sign", "measurement_sign"),
    "measurement_value": ("measurement value", "measurement_value"),
    "measurement_units": ("measurement units", "measurement_units", "measurement unit", "measurement_unit"),
    "laboratory_method": (
        "laboratory typing method",
        "laboratory_typing_method",
        "testing method",
    ),
    "testing_standard": ("testing standard", "testing_standard"),
    "testing_standard_year": ("testing standard year", "testing_standard_year"),
    "source_detail": ("source", "publication", "pubmed"),
}

PREDICTED_MARKERS = re.compile(r"\b(?:predict(?:ed|ion)?|computational|machine learning|model)\b", re.I)


def normalize_header(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower().replace("_", " "))


def resolve_columns(columns: Iterable[str]) -> dict[str, str | None]:
    normalized = {normalize_header(column): column for column in columns}
    resolved: dict[str, str | None] = {}
    for semantic_name, candidates in COLUMN_CANDIDATES.items():
        resolved[semantic_name] = next(
            (normalized[normalize_header(candidate)] for candidate in candidates if normalize_header(candidate) in normalized),
            None,
        )
    required = ("genome_id", "genome_name", "antibiotic", "phenotype")
    missing = [name for name in required if resolved[name] is None]
    if missing:
        raise ValueError(f"BV-BRC table is missing required fields {missing}; columns={list(columns)}")
    return resolved


def normalize_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.replace(r"\s+", " ", regex=True)


def normalize_antibiotic(value: str) -> str:
    key = re.sub(r"\s+", " ", str(value).strip().lower())
    return cfg.ANTIBIOTIC_ALIASES.get(key, key)


def species_mask(df: pd.DataFrame, columns: Mapping[str, str | None], species: str) -> pd.Series:
    species_column = columns.get("species")
    if species_column:
        observed = normalize_text(df[species_column]).str.lower()
        return observed.eq(species.lower())
    names = normalize_text(df[columns["genome_name"]])
    pattern = rf"^{re.escape(species)}(?:\s|$)"
    return names.str.match(pattern, case=False, na=False)


def is_lab_observation(df: pd.DataFrame, columns: Mapping[str, str | None]) -> pd.Series:
    evidence_column = columns.get("evidence")
    if evidence_column is None:
        # The official PATRIC_genome_AMR table contract identifies these rows as
        # laboratory phenotype observations. Absence of a field is recorded in
        # the manifest rather than silently invented per row.
        return pd.Series(True, index=df.index)
    evidence = normalize_text(df[evidence_column])
    return ~evidence.str.contains(PREDICTED_MARKERS, na=False)


def extract_observations(
    raw: pd.DataFrame,
    columns: Mapping[str, str | None],
    targets: Iterable[tuple[str, str]],
) -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    antibiotic_values = normalize_text(raw[columns["antibiotic"]]).map(normalize_antibiotic)
    phenotype_values = normalize_text(raw[columns["phenotype"]]).str.lower()
    lab_mask = is_lab_observation(raw, columns)

    for species, antibiotic in targets:
        target_antibiotic = normalize_antibiotic(antibiotic)
        selected = raw[
            species_mask(raw, columns, species)
            & antibiotic_values.eq(target_antibiotic)
            & phenotype_values.isin(("resistant", "susceptible"))
            & lab_mask
        ].copy()
        if selected.empty:
            continue

        out = pd.DataFrame(index=selected.index)
        out["genome_id"] = normalize_text(selected[columns["genome_id"]])
        out["genome_name"] = normalize_text(selected[columns["genome_name"]])
        out["species"] = species
        out["antibiotic"] = target_antibiotic
        out["phenotype"] = phenotype_values.loc[selected.index]
        out["source"] = "bvbrc"
        out["label_type"] = "lab_measured"
        for field in (
            "evidence",
            "measurement",
            "measurement_sign",
            "measurement_value",
            "measurement_units",
            "laboratory_method",
            "testing_standard",
            "testing_standard_year",
            "source_detail",
        ):
            source_column = columns.get(field)
            out[field] = normalize_text(selected[source_column]) if source_column else ""
        records.append(out.reset_index(drop=True))

    if not records:
        return pd.DataFrame()
    observations = pd.concat(records, ignore_index=True)
    return observations[observations["genome_id"].ne("")].drop_duplicates().reset_index(drop=True)


def separate_conflicts(observations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if observations.empty:
        return observations.copy(), observations.copy()
    keys = ["genome_id", "species", "antibiotic"]
    phenotype_counts = observations.groupby(keys)["phenotype"].transform("nunique")
    conflicts = observations[phenotype_counts.gt(1)].copy()
    clean = observations[phenotype_counts.eq(1)].copy()
    return clean, conflicts


def collapse_labels(clean_observations: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "genome_id",
        "species",
        "antibiotic",
        "phenotype",
        "source",
        "label_type",
        "observation_count",
    ]
    if clean_observations.empty:
        return pd.DataFrame(columns=columns)
    keys = ["genome_id", "species", "antibiotic", "phenotype", "source", "label_type"]
    return (
        clean_observations.groupby(keys, as_index=False)
        .size()
        .rename(columns={"size": "observation_count"})[columns]
        .sort_values(["species", "antibiotic", "genome_id"])
        .reset_index(drop=True)
    )


def evaluate_gate(labels: pd.DataFrame, targets: Iterable[tuple[str, str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for species, antibiotic in targets:
        canonical_antibiotic = normalize_antibiotic(antibiotic)
        pair = labels[
            labels["species"].eq(species) & labels["antibiotic"].eq(canonical_antibiotic)
        ]
        counts = pair["phenotype"].value_counts()
        resistant = int(counts.get("resistant", 0))
        susceptible = int(counts.get("susceptible", 0))
        total = resistant + susceptible
        minority_fraction = min(resistant, susceptible) / total if total else 0.0
        problems: list[str] = []
        if total < cfg.GATE["min_total"]:
            problems.append(f"total {total} < {cfg.GATE['min_total']}")
        if resistant < cfg.GATE["min_per_class"] or susceptible < cfg.GATE["min_per_class"]:
            problems.append(
                f"class count R={resistant}, S={susceptible} < {cfg.GATE['min_per_class']}"
            )
        if minority_fraction < cfg.GATE["min_minority_frac"]:
            problems.append(
                f"minority fraction {minority_fraction:.3f} < {cfg.GATE['min_minority_frac']:.3f}"
            )
        rows.append(
            {
                "species": species,
                "antibiotic": canonical_antibiotic,
                "total": total,
                "resistant": resistant,
                "susceptible": susceptible,
                "minority_fraction": minority_fraction,
                "preliminary_gate_pass": not problems,
                "problems": "; ".join(problems),
                "note": "Genetic-group viability is not evaluated yet.",
            }
        )
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_table(url: str, destination: Path, refresh: bool) -> None:
    if destination.exists() and destination.stat().st_size and not refresh:
        return
    if shutil.which("curl") is None:
        raise RuntimeError("curl is required to retrieve the official FTPS dataset")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    command = [
        "curl",
        "--fail",
        "--location",
        "--silent",
        "--show-error",
        "--connect-timeout",
        "30",
        "--max-time",
        "900",
        "--ssl-reqd",
        "--user",
        "anonymous:guest",
        "--output",
        str(partial),
        url,
    ]
    try:
        subprocess.run(command, check=True)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    if partial.stat().st_size == 0:
        partial.unlink(missing_ok=True)
        raise RuntimeError("BV-BRC returned an empty phenotype table")
    partial.replace(destination)


def download_api_table(destination: Path, refresh: bool) -> None:
    """Fetch only configured species/drug records through the HTTPS Data API.

    This is a practical fallback when the large FTPS release is unreachable.
    The API response is converted to the same TSV shape consumed by ``audit``.
    """

    if destination.exists() and destination.stat().st_size and not refresh:
        return
    if shutil.which("curl") is None:
        raise RuntimeError("curl is required to retrieve the BV-BRC HTTPS API data")
    records: list[dict[str, object]] = []
    for species, antibiotic in cfg.CANDIDATE_TARGETS:
        taxon_id = cfg.TAXON_IDS.get(species)
        if taxon_id is None:
            raise ValueError(f"No BV-BRC taxon ID configured for {species}")
        query = (
            f"eq(taxon_id,{taxon_id})&eq(antibiotic,{quote(antibiotic)})&limit(100000)"
        )
        url = f"{AMR_API_URL}?{query}"
        command = [
            "curl", "--fail", "--location", "--silent", "--show-error",
            "--connect-timeout", "30", "--max-time", "900", url,
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        payload = json.loads(result.stdout)
        if not isinstance(payload, list):
            raise RuntimeError(f"Unexpected BV-BRC API response for {species}/{antibiotic}")
        records.extend(record for record in payload if isinstance(record, dict))
    if not records:
        raise RuntimeError("BV-BRC API returned no records for the configured candidates")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    pd.DataFrame(records).to_csv(partial, sep="\t", index=False)
    partial.replace(destination)


def atomic_csv(frame: pd.DataFrame, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    frame.to_csv(partial, index=False)
    partial.replace(destination)


def audit(source_table: Path, out_dir: Path) -> dict[str, object]:
    raw = pd.read_csv(source_table, sep="\t", dtype=str, low_memory=False)
    columns = resolve_columns(raw.columns)
    observations = extract_observations(raw, columns, cfg.CANDIDATE_TARGETS)
    clean_observations, conflicts = separate_conflicts(observations)
    labels = collapse_labels(clean_observations)
    gate_report = evaluate_gate(labels, cfg.CANDIDATE_TARGETS)

    viable_pairs = {
        (row["species"], row["antibiotic"])
        for _, row in gate_report.iterrows()
        if bool(row["preliminary_gate_pass"])
    }
    viable_labels = labels[
        labels.apply(lambda row: (row["species"], row["antibiotic"]) in viable_pairs, axis=1)
    ].reset_index(drop=True)

    atomic_csv(observations, out_dir / "normalized_observations.csv")
    atomic_csv(conflicts, out_dir / "conflicting_observations.csv")
    atomic_csv(labels, out_dir / "all_labels.csv")
    atomic_csv(viable_labels, out_dir / "labels.csv")
    atomic_csv(gate_report, out_dir / "gate_report.csv")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "BV-BRC PATRIC_genome_AMR",
        "source_file": str(source_table),
        "source_sha256": sha256(source_table),
        "source_bytes": source_table.stat().st_size,
        "resolved_columns": columns,
        "lab_evidence_basis": (
            f"source field {columns['evidence']} filtered for predicted markers"
            if columns.get("evidence")
            else "official PATRIC_genome_AMR table contract; no evidence column present"
        ),
        "candidate_targets": list(cfg.CANDIDATE_TARGETS),
        "gate": cfg.GATE,
        "normalized_observations": len(observations),
        "conflicting_observations": len(conflicts),
        "labels_after_conflict_removal": len(labels),
        "viable_labels_forwarded": len(viable_labels),
    }
    manifest_path = out_dir / "provenance_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {"manifest": manifest, "gate_report": gate_report}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        help="Audit an existing PATRIC_genome_AMR.txt instead of downloading it.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--source-url", default=AMR_TABLE_URL)
    parser.add_argument(
        "--api",
        action="store_true",
        help="Use the HTTPS BV-BRC Data API for configured candidates instead of the FTPS release.",
    )
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--labels-only",
        action="store_true",
        help="Compatibility flag; this command always stops after the label audit.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source_table = args.input.resolve() if args.input else out_dir / (
        "PATRIC_genome_AMR_api.tsv" if args.api else "PATRIC_genome_AMR.txt"
    )
    if args.input is None:
        if args.api:
            print(f"[fetch] {AMR_API_URL}")
            download_api_table(source_table, args.refresh)
        else:
            print(f"[fetch] {args.source_url}")
            download_table(args.source_url, source_table, args.refresh)
    print(f"[audit] {source_table}")
    result = audit(source_table, out_dir)
    print(result["gate_report"].to_string(index=False))
    print(f"[write] receipts: {out_dir}")
    print("[halt] Human Gate #1: review provenance, conflicts, and preliminary viability.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

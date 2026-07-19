from __future__ import annotations

import argparse
import textwrap
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
LABELS_PATH = ROOT / "data" / "processed" / "bvbrc_ecoli_labels.csv"
COHORT_LABELS_PATH = ROOT / "data" / "processed" / "ecoli_cohort_labels.csv"
METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
SKIPPED_FASTA_PATH = ROOT / "reports" / "metrics" / "ecoli_skipped_fasta.csv"
FASTA_DIR = ROOT / "data" / "raw" / "fasta"
API_BASE = "https://www.bv-brc.org/api"

DEFAULT_ANTIBIOTICS = [
    "ampicillin",
    "amoxicillin",
    "cefotaxime",
    "gentamicin",
    "chloramphenicol",
]


def choose_antibiotics(labels: pd.DataFrame, requested: str | None, top_n: int) -> list[str]:
    if requested:
        return [value.strip().lower() for value in requested.split(",") if value.strip()]

    counts = labels.pivot_table(
        index="antibiotic",
        columns="label",
        values="genome_id",
        aggfunc="count",
        fill_value=0,
    )
    for label in ["resistant", "susceptible"]:
        if label not in counts.columns:
            counts[label] = 0
    counts["minority"] = counts[["resistant", "susceptible"]].min(axis=1)
    counts["total"] = counts[["resistant", "susceptible"]].sum(axis=1)
    usable = counts[counts["minority"] >= 10].sort_values(["minority", "total"], ascending=False)
    if usable.empty:
        return DEFAULT_ANTIBIOTICS
    return list(usable.head(top_n).index)


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def fetch_genome_metadata(genome_ids: list[str]) -> pd.DataFrame:
    genome_ids = [str(value) for value in genome_ids]
    if not genome_ids:
        return pd.DataFrame()
    fields = ",".join(
        [
            "genome_id",
            "genome_name",
            "species",
            "genome_quality",
            "genome_status",
            "genome_length",
            "contigs",
            "checkm_completeness",
            "checkm_contamination",
            "cgmlst_hc100",
            "cgmlst_hc50",
            "cgmlst_hc20",
            "isolation_country",
            "host_common_name",
            "collection_year",
        ]
    )
    rows: list[dict] = []
    session = requests.Session()
    for ids in chunked(genome_ids, 100):
        encoded = ",".join(quote(value) for value in ids)
        url = (
            f"{API_BASE}/genome/?in(genome_id,({encoded}))"
            f"&select({fields})&limit({len(ids)})&http_accept=application/json"
        )
        response = session.get(url, timeout=60)
        response.raise_for_status()
        rows.extend(response.json())
    return pd.DataFrame(rows)


def filter_quality_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    if metadata.empty:
        return metadata

    frame = metadata.copy()
    frame["genome_id"] = frame["genome_id"].astype(str)
    contamination = pd.to_numeric(frame.get("checkm_contamination"), errors="coerce")
    completeness = pd.to_numeric(frame.get("checkm_completeness"), errors="coerce")
    quality = frame.get("genome_quality", "").fillna("").astype(str).str.lower()

    keep = quality.eq("good")
    keep &= contamination.isna() | contamination.le(5)
    keep &= completeness.isna() | completeness.ge(95)
    return frame[keep].copy()


def choose_balanced_genomes(labels: pd.DataFrame, max_genomes: int) -> list[str]:
    labels = labels.copy()
    labels["genome_id"] = labels["genome_id"].astype(str)
    counts = labels.groupby("genome_id")["antibiotic"].nunique().sort_values(ascending=False)
    candidates = list(counts.index)
    selected: list[str] = []
    per_label_seen: dict[tuple[str, str], int] = {}

    # Greedy pass: prefer genomes that contribute rare antibiotic/label combinations.
    scored = []
    pair_counts = labels.groupby(["antibiotic", "label"]).size().to_dict()
    for genome_id, group in labels.groupby("genome_id"):
        score = 0.0
        for row in group.itertuples():
            score += 1.0 / max(pair_counts[(row.antibiotic, row.label)], 1)
        scored.append((score, genome_id))
    for _, genome_id in sorted(scored, reverse=True):
        genome_rows = labels[labels["genome_id"] == genome_id]
        improves_balance = False
        for row in genome_rows.itertuples():
            key = (row.antibiotic, row.label)
            if per_label_seen.get(key, 0) < max(5, max_genomes // 20):
                improves_balance = True
                break
        if improves_balance or genome_id in candidates[:max_genomes]:
            selected.append(genome_id)
            for row in genome_rows.itertuples():
                key = (row.antibiotic, row.label)
                per_label_seen[key] = per_label_seen.get(key, 0) + 1
        if len(selected) >= max_genomes:
            break
    return selected


def existing_fasta_ids() -> set[str]:
    if not FASTA_DIR.exists():
        return set()
    return {path.stem for path in FASTA_DIR.glob("*.fna") if path.stat().st_size > 0}


def download_fasta(genome_id: str, output_path: Path) -> None:
    if output_path.exists() and output_path.stat().st_size > 0:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    url = (
        f"{API_BASE}/genome_sequence/?eq(genome_id,{quote(genome_id)})"
        "&select(sequence_id,description,sequence)&limit(25000)"
        "&http_accept=application/json"
    )
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    records = response.json()
    if not records:
        raise RuntimeError(f"No sequence records found for genome {genome_id}")

    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            sequence_id = record.get("sequence_id") or f"{genome_id}_unknown"
            description = record.get("description") or ""
            sequence = (record.get("sequence") or "").upper()
            handle.write(f">{sequence_id} {description}".rstrip() + "\n")
            handle.write("\n".join(textwrap.wrap(sequence, width=80)) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-genomes", type=int, default=200)
    parser.add_argument("--target-genomes", type=int, default=50)
    parser.add_argument("--antibiotics", help="Comma-separated antibiotics. Defaults to auto-selected balanced drugs.")
    parser.add_argument("--top-antibiotics", type=int, default=5)
    parser.add_argument("--skip-fasta", action="store_true")
    parser.add_argument(
        "--use-existing-fasta",
        action="store_true",
        help="Build the cohort only from FASTA files already present in data/raw/fasta.",
    )
    args = parser.parse_args()

    labels = pd.read_csv(LABELS_PATH)
    labels["genome_id"] = labels["genome_id"].astype(str)
    if args.use_existing_fasta:
        fasta_ids = existing_fasta_ids()
        print(f"Existing non-empty FASTA files: {len(fasta_ids)}", flush=True)
        labels = labels[labels["genome_id"].isin(fasta_ids)].copy()
        if labels.empty:
            raise SystemExit("No BV-BRC labels match existing FASTA files in data/raw/fasta.")

    selected_antibiotics = choose_antibiotics(labels, args.antibiotics, args.top_antibiotics)
    print(f"Selected antibiotics: {', '.join(selected_antibiotics)}", flush=True)
    labels = labels[labels["antibiotic"].isin(selected_antibiotics)].copy()
    labels = labels.drop_duplicates(["genome_id", "antibiotic", "label"])

    candidate_ids = choose_balanced_genomes(labels, args.max_genomes)
    print(f"Candidate genomes selected before metadata filter: {len(candidate_ids)}", flush=True)
    candidate_metadata = fetch_genome_metadata(candidate_ids)
    print(f"Fetched metadata rows: {len(candidate_metadata)}", flush=True)
    quality_metadata = filter_quality_metadata(candidate_metadata)
    quality_ids = set(quality_metadata["genome_id"].astype(str)) if not quality_metadata.empty else set()
    selected_ids = [genome_id for genome_id in candidate_ids if genome_id in quality_ids]

    available_ids = selected_ids[: args.target_genomes]
    skipped_rows: list[dict[str, str]] = []
    if args.use_existing_fasta:
        fasta_ids = existing_fasta_ids()
        available_ids = [genome_id for genome_id in selected_ids if genome_id in fasta_ids][: args.target_genomes]
        missing_after_quality = [genome_id for genome_id in selected_ids if genome_id not in fasta_ids]
        skipped_rows.extend(
            {"genome_id": genome_id, "reason": "quality-passing metadata but local FASTA missing"}
            for genome_id in missing_after_quality
        )
    elif not args.skip_fasta:
        available_ids = []
        for index, genome_id in enumerate(selected_ids, start=1):
            print(f"Downloading FASTA {index}/{len(selected_ids)}: {genome_id}", flush=True)
            try:
                download_fasta(genome_id, FASTA_DIR / f"{genome_id}.fna")
                available_ids.append(genome_id)
            except Exception as exc:
                print(f"Skipping {genome_id}: {exc}", flush=True)
                skipped_rows.append({"genome_id": genome_id, "reason": str(exc)})
            if len(available_ids) >= args.target_genomes:
                break

    cohort_labels = labels[labels["genome_id"].isin(available_ids)].copy()
    metadata = quality_metadata[quality_metadata["genome_id"].isin(available_ids)].copy()

    COHORT_LABELS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SKIPPED_FASTA_PATH.parent.mkdir(parents=True, exist_ok=True)
    cohort_labels.to_csv(COHORT_LABELS_PATH, index=False)
    metadata.to_csv(METADATA_PATH, index=False)
    pd.DataFrame(skipped_rows).to_csv(SKIPPED_FASTA_PATH, index=False)

    print(f"Candidate genomes: {len(candidate_ids)}", flush=True)
    print(f"Quality-filtered genomes: {len(selected_ids)}", flush=True)
    print(f"Genomes with FASTA: {len(available_ids)}", flush=True)
    print(f"Cohort labels: {len(cohort_labels)}", flush=True)
    print("Cohort label counts:", flush=True)
    print(
        cohort_labels.pivot_table(
            index="antibiotic",
            columns="label",
            values="genome_id",
            aggfunc="count",
            fill_value=0,
        )
    )
    print(f"Wrote {COHORT_LABELS_PATH}", flush=True)
    print(f"Wrote {METADATA_PATH}", flush=True)


if __name__ == "__main__":
    main()

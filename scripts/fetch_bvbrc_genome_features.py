from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
OUT_PATH = ROOT / "data" / "interim" / "bvbrc_genome_features.csv"
API_BASE = "https://www.bv-brc.org/api"
FIELDS = [
    "genome_id",
    "patric_id",
    "refseq_locus_tag",
    "gene",
    "product",
    "feature_type",
    "annotation",
]


def fetch_features_for_genome(session: requests.Session, genome_id: str, limit: int) -> list[dict]:
    fields = ",".join(FIELDS)
    url = (
        f"{API_BASE}/genome_feature/?eq(genome_id,{quote(genome_id)})"
        f"&select({fields})&limit({limit})&http_accept=application/json"
    )
    response = session.get(url, timeout=120)
    response.raise_for_status()
    rows = response.json()
    if not isinstance(rows, list):
        raise RuntimeError(f"Unexpected BV-BRC response for genome {genome_id}: {type(rows).__name__}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch BV-BRC genome_feature annotations for the current E. coli cohort."
    )
    parser.add_argument("--metadata", type=Path, default=METADATA_PATH)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    parser.add_argument("--limit-per-genome", type=int, default=25000)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    genome_ids = sorted({str(value) for value in metadata["genome_id"].dropna()})
    existing = pd.DataFrame()
    done_ids: set[str] = set()
    if args.out.exists() and not args.force:
        existing = pd.read_csv(args.out, dtype=str)
        if "genome_id" in existing.columns:
            done_ids = set(existing["genome_id"].dropna().astype(str))
        print(f"Existing feature rows: {len(existing)} for {len(done_ids)} genomes", flush=True)

    rows: list[dict] = []
    session = requests.Session()
    pending_ids = [genome_id for genome_id in genome_ids if genome_id not in done_ids]
    for index, genome_id in enumerate(pending_ids, start=1):
        print(f"Fetching BV-BRC features {index}/{len(pending_ids)}: {genome_id}", flush=True)
        try:
            fetched = fetch_features_for_genome(session, genome_id, args.limit_per_genome)
        except Exception as exc:
            print(f"Skipping feature fetch for {genome_id}: {exc}", flush=True)
            rows.append({"genome_id": genome_id, "fetch_error": str(exc)})
            continue
        if fetched:
            rows.extend(fetched)
        else:
            rows.append({"genome_id": genome_id, "fetch_error": "no_features_returned"})

    new_frame = pd.DataFrame(rows)
    if not existing.empty and not new_frame.empty:
        frame = pd.concat([existing, new_frame], ignore_index=True, sort=False)
    elif not existing.empty:
        frame = existing
    else:
        frame = new_frame

    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(frame)} rows", flush=True)


if __name__ == "__main__":
    main()

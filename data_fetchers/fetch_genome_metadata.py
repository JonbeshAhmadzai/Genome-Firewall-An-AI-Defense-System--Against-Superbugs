#!/usr/bin/env python3
"""Fetch quality and lineage metadata for an explicit BV-BRC genome cohort."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


API_BASE = "https://www.bv-brc.org/api/genome/"
FIELDS = (
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
)


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-list", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("data/processed/cohort100_genome_metadata.csv"))
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()
    genome_ids = list(
        dict.fromkeys(
            line.strip()
            for line in args.genome_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    )
    fields = ",".join(FIELDS)
    rows: list[dict] = []
    session = requests.Session()
    for batch in chunks(genome_ids, args.batch_size):
        encoded = ",".join(quote(value) for value in batch)
        url = f"{API_BASE}?in(genome_id,({encoded}))&select({fields})&limit({len(batch)})&http_accept=application/json"
        response = session.get(url, timeout=120)
        response.raise_for_status()
        payload = response.json()
        if not isinstance(payload, list):
            raise RuntimeError(f"Unexpected BV-BRC metadata response: {type(payload).__name__}")
        rows.extend(payload)
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["genome_id"] = frame["genome_id"].astype(str)
        frame = frame.drop_duplicates("genome_id").set_index("genome_id").reindex(genome_ids).reset_index()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(frame)} rows for {len(genome_ids)} requested genomes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

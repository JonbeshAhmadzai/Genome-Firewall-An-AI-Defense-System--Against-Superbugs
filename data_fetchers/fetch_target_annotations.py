#!/usr/bin/env python3
"""Fetch resumable BV-BRC general annotations for molecular-target checks."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


API_BASE = "https://www.bv-brc.org/api/genome_feature/"
FIELDS = ("genome_id", "patric_id", "refseq_locus_tag", "gene", "product", "feature_type", "annotation")


def fetch_one(session: requests.Session, genome_id: str, limit: int) -> pd.DataFrame:
    fields = ",".join(FIELDS)
    url = (
        f"{API_BASE}?eq(genome_id,{quote(genome_id)})"
        f"&select({fields})&limit({limit})&http_accept=application/json"
    )
    response = session.get(url, timeout=120)
    response.raise_for_status()
    rows = response.json()
    if not isinstance(rows, list) or not rows:
        raise RuntimeError("no general genome annotations returned")
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-list", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("data/raw/bvbrc/cohort200/target_annotations"))
    parser.add_argument("--limit-per-genome", type=int, default=25000)
    args = parser.parse_args()
    genome_ids = list(
        dict.fromkeys(
            line.strip()
            for line in args.genome_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    session = requests.Session()
    for index, genome_id in enumerate(genome_ids, 1):
        output = args.out_dir / f"{genome_id}.tsv.gz"
        if output.exists() and output.stat().st_size > 0:
            status = "existing"
            error = ""
        else:
            partial = output.with_suffix(output.suffix + ".part")
            try:
                frame = fetch_one(session, genome_id, args.limit_per_genome)
                frame.to_csv(partial, sep="\t", index=False, compression="gzip")
                partial.replace(output)
                status = "downloaded"
                error = ""
            except Exception as exc:
                partial.unlink(missing_ok=True)
                status = "failed"
                error = str(exc)
        rows.append({"genome_id": genome_id, "status": status, "output": str(output), "error": error})
        print(f"[{index}/{len(genome_ids)}] {genome_id}: {status}", flush=True)
        (args.out_dir / "fetch_manifest.json").write_text(
            json.dumps(
                {
                    "updated_at_utc": datetime.now(timezone.utc).isoformat(),
                    "requested": len(genome_ids),
                    "completed": sum(row["status"] in {"downloaded", "existing"} for row in rows),
                    "failed": sum(row["status"] == "failed" for row in rows),
                    "rows": rows,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return 0 if all(row["status"] != "failed" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())

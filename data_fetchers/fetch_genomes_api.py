#!/usr/bin/env python3
"""Download selected BV-BRC assemblies through the HTTPS genome-sequence API."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote


API_URL = "https://www.bv-brc.org/api/genome_sequence/"
DEFAULT_LIST = Path("data/raw/bvbrc/selected/genome_list.txt")
DEFAULT_OUT = Path("data/raw/bvbrc/selected/genomes")


def fetch_records(genome_ids: list[str]) -> list[dict[str, object]]:
    joined = ",".join(quote(genome_id) for genome_id in genome_ids)
    query = f"in(genome_id,({joined}))&limit(100000)"
    command = [
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--connect-timeout", "30", "--max-time", "900",
        f"{API_URL}?{query}",
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    records = json.loads(result.stdout)
    if not isinstance(records, list) or not records:
        raise ValueError(f"No contigs returned for genomes {genome_ids}")
    return [record for record in records if isinstance(record, dict) and record.get("sequence")]


def write_genome_records(genome_id: str, records: list[dict[str, object]], output: Path) -> dict[str, object]:
    records = [record for record in records if str(record.get("genome_id")) == genome_id]
    if not records:
        raise ValueError(f"No contigs returned for genome {genome_id}")
    records.sort(key=lambda record: str(record.get("sequence_id", "")))
    total_bases = 0
    partial = output.with_suffix(output.suffix + ".part")
    partial.parent.mkdir(parents=True, exist_ok=True)
    import gzip

    fasta_chunks: list[str] = []
    for record in records:
        sequence_id = str(record.get("sequence_id") or record.get("accession"))
        sequence = str(record["sequence"]).upper()
        total_bases += len(sequence)
        fasta_chunks.append(f">{sequence_id} genome={genome_id}\n")
        fasta_chunks.extend(sequence[start : start + 80] + "\n" for start in range(0, len(sequence), 80))
    with gzip.open(partial, "wb") as stream:
        stream.write("".join(fasta_chunks).encode("ascii"))
    partial.replace(output)
    return {
        "genome_id": genome_id,
        "contigs": len(records),
        "total_bases": total_bases,
        "bytes": output.stat().st_size,
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "status": "downloaded",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-list", type=Path, default=DEFAULT_LIST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-genomes", type=int, default=20)
    parser.add_argument("--min-free-gb", type=float, default=2.0)
    args = parser.parse_args()
    ids = [line.strip() for line in args.genome_list.read_text().splitlines() if line.strip()]
    # Preserve the selector's intentional balanced ordering.
    ids = list(dict.fromkeys(ids))[: args.max_genomes]
    if not ids:
        raise SystemExit(f"No genome IDs found in {args.genome_list}")
    free_gb = shutil.disk_usage(args.out_dir.parent if args.out_dir.parent.exists() else Path.cwd()).free / (1024**3)
    if free_gb < args.min_free_gb:
        raise SystemExit(f"Only {free_gb:.2f} GB free; refusing to download with less than {args.min_free_gb:.2f} GB reserved")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = args.out_dir / "download_receipt.csv"
    existing = {}
    if receipt_path.exists():
        with receipt_path.open(newline="") as stream:
            existing = {row["genome_id"]: row for row in csv.DictReader(stream)}
    rows = []
    pending: list[tuple[int, str, Path]] = []
    for index, genome_id in enumerate(ids, 1):
        destination = args.out_dir / f"{genome_id}.fna.gz"
        if destination.exists() and destination.stat().st_size > 0:
            rows.append(existing.get(genome_id, {"genome_id": genome_id, "status": "existing", "bytes": destination.stat().st_size}))
            continue
        pending.append((index, genome_id, destination))

    batch_size = 10
    for start in range(0, len(pending), batch_size):
        batch = pending[start : start + batch_size]
        try:
            records = fetch_records([genome_id for _, genome_id, _ in batch])
            by_genome = {}
            for record in records:
                by_genome.setdefault(str(record.get("genome_id")), []).append(record)
            for index, genome_id, destination in batch:
                try:
                    receipt = write_genome_records(genome_id, by_genome.get(genome_id, []), destination)
                except Exception as error:
                    receipt = {"genome_id": genome_id, "status": "failed", "error": str(error)}
                    destination.unlink(missing_ok=True)
                rows.append(receipt)
                print(f"[{index}/{len(ids)}] {genome_id}: {receipt['status']}")
        except Exception as error:
            for index, genome_id, destination in batch:
                rows.append({"genome_id": genome_id, "status": "failed", "error": str(error)})
                destination.unlink(missing_ok=True)
                print(f"[{index}/{len(ids)}] {genome_id}: failed")
        fields = sorted({key for row in rows for key in row})
        with receipt_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    fields = sorted({key for row in rows for key in row})
    with receipt_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "api": API_URL,
        "requested": len(ids),
        "downloaded": sum(row.get("status") == "downloaded" for row in rows),
        "failed": sum(row.get("status") == "failed" for row in rows),
        "output": str(args.out_dir),
    }
    (args.out_dir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0 if manifest["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

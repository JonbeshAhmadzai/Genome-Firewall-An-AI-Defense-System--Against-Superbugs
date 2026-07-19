#!/usr/bin/env python3
"""Run AMRFinderPlus over selected compressed assemblies, resumably."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.genome_reader.run_amrfinder import run_amrfinder


def genome_id(path: Path) -> str:
    return path.name[:-len(".fna.gz")] if path.name.endswith(".fna.gz") else path.stem


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", type=Path, default=Path("data/raw/bvbrc/selected/genomes"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed/amrfinder"))
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--genome-list", type=Path, default=None)
    parser.add_argument("--executable", default="amrfinder")
    parser.add_argument("--organism", default="Escherichia")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--plus", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    paths = sorted(args.genome_dir.glob("*.fna.gz"))
    if args.genome_list is not None:
        requested = {
            line.strip()
            for line in args.genome_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }
        paths = [path for path in paths if genome_id(path) in requested]
    if args.max_genomes is not None:
        paths = paths[: args.max_genomes]
    if not paths:
        raise SystemExit(f"No compressed FASTAs found in {args.genome_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, path in enumerate(paths, 1):
        gid = genome_id(path)
        output = args.out_dir / f"{gid}.tsv"
        if output.exists() and output.stat().st_size > 0:
            row = {"genome_id": gid, "status": "existing", "output": str(output)}
        else:
            try:
                run_amrfinder(
                    path,
                    output,
                    executable=args.executable,
                    organism=args.organism,
                    plus=args.plus,
                    threads=args.threads,
                )
                row = {"genome_id": gid, "status": "annotated", "output": str(output)}
            except Exception as error:
                output.unlink(missing_ok=True)
                row = {"genome_id": gid, "status": "failed", "error": str(error)}
        rows.append(row)
        print(f"[{index}/{len(paths)}] {gid}: {row['status']}", flush=True)
    receipt = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested": len(paths),
        "annotated": sum(row["status"] in {"annotated", "existing"} for row in rows),
        "failed": sum(row["status"] == "failed" for row in rows),
        "rows": rows,
    }
    (args.out_dir / "annotation_manifest.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return 0 if receipt["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

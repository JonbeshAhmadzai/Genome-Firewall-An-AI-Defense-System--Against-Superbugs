#!/usr/bin/env python3
"""Run AMRFinderPlus over selected compressed assemblies, resumably."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.genome_reader.run_amrfinder import run_amrfinder
from targets_config import enabled_species, pathogen_spec


def genome_id(path: Path) -> str:
    return path.name[:-len(".fna.gz")] if path.name.endswith(".fna.gz") else path.stem


def annotate_one(
    path: Path,
    *,
    out_dir: Path,
    executable: str,
    organism: str | None,
    plus: bool,
    threads: int,
) -> dict[str, str]:
    gid = genome_id(path)
    output = out_dir / f"{gid}.tsv"
    if output.exists() and output.stat().st_size > 0:
        return {"genome_id": gid, "status": "existing", "output": str(output)}
    try:
        run_amrfinder(
            path,
            output,
            executable=executable,
            organism=organism,
            plus=plus,
            threads=threads,
        )
        return {"genome_id": gid, "status": "annotated", "output": str(output)}
    except Exception as error:
        output.unlink(missing_ok=True)
        return {"genome_id": gid, "status": "failed", "error": str(error)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", type=Path, default=Path("data/raw/bvbrc/selected/genomes"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed/amrfinder"))
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--genome-list", type=Path, default=None)
    parser.add_argument("--executable", default="amrfinder")
    parser.add_argument(
        "--species",
        default=enabled_species(kind="bacterium")[0],
        help="Configured bacterial species; its AMRFinderPlus organism is selected automatically.",
    )
    parser.add_argument("--organism", default=None, help="Optional AMRFinderPlus override.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--workers", type=int, default=1, help="Concurrent AMRFinder jobs; use a bounded value such as 4.")
    parser.add_argument("--plus", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    if args.threads < 1 or args.workers < 1:
        raise SystemExit("threads and workers must be positive")
    species_config = pathogen_spec(args.species)
    if species_config.get("reader") != "amrfinderplus":
        raise SystemExit(
            f"Target annotation backend {species_config.get('reader')!r} is not implemented for {args.species}."
        )
    organism = args.organism or species_config.get("amrfinder_organism")
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
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                annotate_one,
                path,
                out_dir=args.out_dir,
                executable=args.executable,
                organism=organism,
                plus=args.plus,
                threads=args.threads,
            ): path
            for path in paths
        }
        for index, future in enumerate(as_completed(futures), 1):
            row = future.result()
            rows.append(row)
            print(f"[{index}/{len(paths)}] {row['genome_id']}: {row['status']}", flush=True)
    rows.sort(key=lambda row: row["genome_id"])
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

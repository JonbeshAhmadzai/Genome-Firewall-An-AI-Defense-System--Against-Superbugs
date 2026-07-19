#!/usr/bin/env python3
"""Create FASTA QC and homology-group receipts for selected assemblies."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.genome_reader.fasta_qc import fasta_stats
from src.genome_reader.grouping import assign_homology_groups


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", type=Path, default=Path("data/raw/bvbrc/selected/genomes"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--genome-list", type=Path, default=None)
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument("--jaccard-threshold", type=float, default=0.90)
    args = parser.parse_args()
    paths = sorted(args.genome_dir.glob("*.fna.gz"))
    if args.genome_list is not None:
        requested = {
            line.strip()
            for line in args.genome_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }
        paths = [
            path
            for path in paths
            if (path.name[:-len(".fna.gz")] if path.name.endswith(".fna.gz") else path.stem) in requested
        ]
    if not paths:
        raise SystemExit(f"No FASTA files found in {args.genome_dir}")
    stats = []
    failed = []
    valid_paths = []
    for path in paths:
        try:
            stats.append(fasta_stats(path))
            valid_paths.append(path)
        except Exception as error:
            failed.append({"path": str(path), "error": str(error)})
    groups = assign_homology_groups(valid_paths, jaccard_threshold=args.jaccard_threshold)
    group_frame = pd.DataFrame({"genome_id": list(groups), "homology_group_id": list(groups.values())})
    if args.metadata is not None and args.metadata.exists():
        metadata = pd.read_csv(args.metadata, dtype=str)
        metadata["genome_id"] = metadata["genome_id"].astype(str)
        if "cgmlst_hc100" in metadata.columns:
            group_frame = group_frame.merge(
                metadata[["genome_id", "cgmlst_hc100"]], on="genome_id", how="left"
            )
            lineage = group_frame["cgmlst_hc100"].fillna("").astype(str).str.strip()
            group_frame.loc[lineage.ne(""), "homology_group_id"] = "cgmlst_hc100:" + lineage[lineage.ne("")]
            group_frame = group_frame.drop(columns=["cgmlst_hc100"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(stats).to_csv(args.out_dir / "genome_qc.csv", index=False)
    group_frame.to_csv(args.out_dir / "genome_groups.csv", index=False)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "jaccard_threshold": args.jaccard_threshold,
        "k": 21,
        "sketch_size": 256,
        "requested": len(paths),
        "valid": len(valid_paths),
        "failed": failed,
        "groups": int(group_frame["homology_group_id"].nunique()),
        "group_source": "cgmlst_hc100_with_sequence_sketch_fallback" if args.metadata else "sequence_jaccard_sketch",
    }
    (args.out_dir / "grouping_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())

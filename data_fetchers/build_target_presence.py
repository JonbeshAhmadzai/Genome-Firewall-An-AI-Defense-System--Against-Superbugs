#!/usr/bin/env python3
"""Build the cohort molecular-target table from fetched BV-BRC annotations."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.genome_reader.target_presence import build_target_presence
from targets_config import MOLECULAR_TARGETS, TARGET_ANNOTATION_MARKERS


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-list", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, default=Path("data/raw/bvbrc/cohort200/target_annotations"))
    parser.add_argument("--out", type=Path, default=Path("data/processed/target_presence.csv"))
    args = parser.parse_args()
    genome_ids = list(
        dict.fromkeys(
            line.strip()
            for line in args.genome_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    )
    annotations = {}
    for genome_id in genome_ids:
        path = args.annotation_dir / f"{genome_id}.tsv.gz"
        if path.exists() and path.stat().st_size > 0:
            annotations[genome_id] = pd.read_csv(path, sep="\t", dtype=str)
    frame = build_target_presence(genome_ids, annotations, TARGET_ANNOTATION_MARKERS, MOLECULAR_TARGETS)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(frame.groupby(["antibiotic", "target_status"]).size().to_string())
    print(f"Wrote {args.out} with {len(frame)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

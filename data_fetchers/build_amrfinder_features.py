#!/usr/bin/env python3
"""Build a sparse AMRFinderPlus feature/evidence pair for a cohort."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.genome_reader.build_features import build_feature_matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinder-dir", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    paths = sorted(args.amrfinder_dir.glob("*.tsv"))
    if not paths:
        raise SystemExit(f"No AMRFinderPlus TSV files found in {args.amrfinder_dir}")
    matrix, evidence = build_feature_matrix(paths)
    args.features.parent.mkdir(parents=True, exist_ok=True)
    args.evidence.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.features, index=False)
    evidence.to_csv(args.evidence, index=False)
    print(f"Wrote {args.features}: {len(matrix)} genomes, {max(len(matrix.columns) - 1, 0)} features")
    print(f"Wrote {args.evidence}: {len(evidence)} evidence rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

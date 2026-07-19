#!/usr/bin/env python3
"""Select a deterministic, storage-safe genome cohort after Human Gate #1.

The full BV-BRC label table contains more genomes than this workstation can
store as assemblies. This selector keeps every ceftriaxone-labelled genome
(the smallest passing target), then adds class-balanced genomes for the other
passing targets until the requested cap is reached.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from targets_config import ACTIVE_TARGETS


DEFAULT_LABELS = Path("data/raw/bvbrc/labels.csv")
DEFAULT_OUT = Path("data/raw/bvbrc/selected")
# The active list comes from targets_config.py. Add/enable a pair there instead
# of editing this selector when a new validated target is ready.
PASSING_TARGETS = tuple(dict.fromkeys(antibiotic for _, antibiotic in ACTIVE_TARGETS))


def stable_key(genome_id: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}:{genome_id}".encode("utf-8")).hexdigest()


def select(labels: pd.DataFrame, *, max_genomes: int, min_per_class: int, seed: int) -> tuple[pd.DataFrame, dict[str, object]]:
    required = {"genome_id", "antibiotic", "phenotype"}
    missing = required - set(labels.columns)
    if missing:
        raise ValueError(f"Labels are missing columns: {sorted(missing)}")
    labels = labels[labels["antibiotic"].isin(PASSING_TARGETS)].copy()
    labels["_key"] = labels["genome_id"].astype(str).map(lambda value: stable_key(value, seed))
    labels = labels.sort_values(["_key", "genome_id"]).drop_duplicates(["genome_id", "antibiotic"])
    selected: set[str] = set()
    selection_order: list[str] = []
    pairs = [
        (antibiotic, phenotype)
        for antibiotic in PASSING_TARGETS
        for phenotype in ("resistant", "susceptible")
    ]
    ids_by_pair = {
        pair: set(
            labels.loc[
                labels["antibiotic"].eq(pair[0]) & labels["phenotype"].eq(pair[1]),
                "genome_id",
            ].astype(str)
        )
        for pair in pairs
    }

    # Round-robin quotas make small cohorts useful for every passing target.
    # Candidates covering several still-unfilled drug/class cells are preferred.
    while len(selected) < max_genomes:
        underfilled = {
            pair
            for pair in pairs
            if len(selected & ids_by_pair[pair]) < min(min_per_class, len(ids_by_pair[pair]))
        }
        if not underfilled:
            break
        progress = False
        for pair in pairs:
            if pair not in underfilled or len(selected) >= max_genomes:
                continue
            candidates = ids_by_pair[pair] - selected
            if not candidates:
                continue
            chosen = min(
                candidates,
                key=lambda genome_id: (
                    -sum(genome_id in ids_by_pair[other] for other in underfilled),
                    stable_key(genome_id, seed),
                    genome_id,
                ),
            )
            selected.add(chosen)
            selection_order.append(chosen)
            progress = True
        if not progress:
            break

    quota_receipt: list[dict[str, object]] = []
    for antibiotic, phenotype in pairs:
        quota_receipt.append(
            {
                "antibiotic": antibiotic,
                "phenotype": phenotype,
                "requested": min_per_class,
                "selected": len(selected & ids_by_pair[(antibiotic, phenotype)]),
            }
        )

    if len(selected) < max_genomes:
        remaining = labels[~labels["genome_id"].astype(str).isin(selected)].drop_duplicates("genome_id")
        additions = remaining.head(max_genomes - len(selected))["genome_id"].astype(str).tolist()
        selected.update(additions)
        selection_order.extend(additions)

    selected_labels = labels[labels["genome_id"].astype(str).isin(selected)].drop(columns="_key").sort_values(["antibiotic", "genome_id"])
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": seed,
        "max_genomes": max_genomes,
        "min_per_class": min_per_class,
        "passing_targets": list(PASSING_TARGETS),
        "selected_genomes": int(selected_labels["genome_id"].nunique()),
        "selected_labels": int(len(selected_labels)),
        "genome_order": selection_order,
        "quota_receipt": quota_receipt,
        "counts": [
            {"antibiotic": antibiotic, "phenotype": phenotype, "n": int(count)}
            for (antibiotic, phenotype), count in selected_labels.groupby(["antibiotic", "phenotype"]).size().items()
        ],
    }
    return selected_labels, manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-genomes", type=int, default=2000)
    parser.add_argument("--min-per-class", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.max_genomes < 1 or args.min_per_class < 1:
        raise SystemExit("max-genomes and min-per-class must be positive")
    labels = pd.read_csv(args.labels, dtype=str)
    selected, manifest = select(labels, max_genomes=args.max_genomes, min_per_class=args.min_per_class, seed=args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out_dir / "selected_labels.csv", index=False)
    pd.Series(manifest["genome_order"], name="genome_id").to_csv(
        args.out_dir / "genome_list.txt", index=False, header=False
    )
    (args.out_dir / "selection_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    print(selected.groupby(["antibiotic", "phenotype"]).size().to_string())
    print(f"Selected genomes: {manifest['selected_genomes']}")
    print(f"Receipts: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Deterministic sketch-based grouping for near-identical assemblies.

This is a lightweight CPU fallback for the homology-group contract. For the
final dataset, the chosen k-mer/sketch settings and threshold must be recorded
and justified, and the resulting groups must be audited before modeling.
"""

from __future__ import annotations

import hashlib
import gzip
from pathlib import Path

from Bio import SeqIO


def _sketch(path: Path, *, k: int = 21, sketch_size: int = 256) -> set[int]:
    stream_factory = gzip.open if path.suffix == ".gz" else Path.open
    with stream_factory(path, "rt") as stream:
        sequence = "".join(str(record.seq).upper() for record in SeqIO.parse(stream, "fasta"))
    if len(sequence) < k:
        return set()
    window_count = len(sequence) - k + 1
    step = max(1, window_count // 10000)
    hashes = {
        int.from_bytes(hashlib.blake2b(sequence[index : index + k].encode(), digest_size=8).digest(), "big")
        for index in range(0, window_count, step)
        if all(base in "ACGT" for base in sequence[index : index + k])
    }
    return set(sorted(hashes)[:sketch_size])


def _jaccard(left: set[int], right: set[int]) -> float:
    if not left or not right:
        return 0.0
    return len(left & right) / len(left | right)


def assign_homology_groups(
    fasta_paths: list[Path],
    *,
    jaccard_threshold: float = 0.90,
    k: int = 21,
    sketch_size: int = 256,
) -> dict[str, str]:
    """Assign deterministic connected-component group IDs to FASTAs."""

    if not 0 < jaccard_threshold <= 1:
        raise ValueError("jaccard_threshold must be in (0, 1]")
    paths = sorted((Path(path) for path in fasta_paths), key=lambda path: path.name)
    def genome_id(path: Path) -> str:
        return path.name[:-len(".fna.gz")] if path.name.endswith(".fna.gz") else path.stem

    sketches = {genome_id(path): _sketch(path, k=k, sketch_size=sketch_size) for path in paths}
    parent = {genome_id: genome_id for genome_id in sketches}

    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    ids = sorted(sketches)
    for index, left in enumerate(ids):
        for right in ids[index + 1 :]:
            if _jaccard(sketches[left], sketches[right]) >= jaccard_threshold:
                union(left, right)

    roots = {genome_id: find(genome_id) for genome_id in ids}
    canonical = {root: f"group-{position:05d}" for position, root in enumerate(sorted(set(roots.values())), 1)}
    return {genome_id: canonical[root] for genome_id, root in roots.items()}

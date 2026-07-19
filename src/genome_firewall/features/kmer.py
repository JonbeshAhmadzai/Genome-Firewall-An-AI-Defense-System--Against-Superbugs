from __future__ import annotations

from collections import Counter


BASES = set("ACGT")


def read_fasta_text(text: str, max_bases: int | None = 500_000) -> str:
    chunks = []
    total = 0
    for raw_line in text.splitlines():
        line = raw_line.strip().upper()
        if not line or line.startswith(">"):
            continue
        line = "".join(base for base in line if base in BASES)
        if max_bases is not None:
            remaining = max_bases - total
            if remaining <= 0:
                break
            line = line[:remaining]
        chunks.append(line)
        total += len(line)
    return "".join(chunks)


def kmer_feature_row(sequence: str, feature_columns: list[str], k: int = 4) -> dict[str, float]:
    counts: Counter[str] = Counter()
    total = 0
    for index in range(0, max(len(sequence) - k + 1, 0)):
        kmer = sequence[index : index + k]
        if set(kmer) <= BASES:
            counts[kmer] += 1
            total += 1
    total = total or 1
    row = {}
    for column in feature_columns:
        kmer = column.removeprefix("kmer_")
        row[column] = counts[kmer] / total
    return row


"""Quality checks and reproducible statistics for assembled FASTA genomes."""

from __future__ import annotations

import hashlib
import gzip
from pathlib import Path

from Bio import SeqIO


DNA_ALPHABET = set("ACGTN")


def fasta_stats(path: Path) -> dict[str, object]:
    """Validate a nucleotide FASTA and return compact QC statistics."""

    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"FASTA is missing or empty: {path}")
    stream_factory = gzip.open if path.suffix == ".gz" else Path.open
    with stream_factory(path, "rt") as stream:
        records = list(SeqIO.parse(stream, "fasta"))
    if not records:
        raise ValueError(f"FASTA contains no records: {path}")
    lengths = []
    ambiguous = 0
    total = 0
    seen_ids: set[str] = set()
    for record in records:
        record_id = record.id.strip()
        if not record_id or record_id in seen_ids:
            raise ValueError(f"FASTA contains an empty or duplicate sequence ID: {record_id!r}")
        seen_ids.add(record_id)
        sequence = str(record.seq).upper()
        invalid = set(sequence) - DNA_ALPHABET
        if invalid:
            raise ValueError(f"FASTA {path} contains invalid nucleotide symbols: {sorted(invalid)}")
        length = len(sequence)
        if length == 0:
            raise ValueError(f"FASTA contains an empty sequence: {record_id}")
        lengths.append(length)
        total += length
        ambiguous += sequence.count("N")
    ordered = sorted(lengths, reverse=True)
    half = total / 2
    running = 0
    n50 = 0
    for length in ordered:
        running += length
        if running >= half:
            n50 = length
            break
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": str(path),
        "sha256": digest,
        "contigs": len(records),
        "total_bases": total,
        "n50": n50,
        "ambiguous_fraction": ambiguous / total if total else 1.0,
    }

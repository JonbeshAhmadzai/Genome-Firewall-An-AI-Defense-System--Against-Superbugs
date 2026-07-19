from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer


ROOT = Path(__file__).resolve().parents[1]
FASTA_DIR = ROOT / "data" / "raw" / "fasta"
METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
FEATURES_PATH = ROOT / "data" / "interim" / "ecoli_kmer_features.csv"
BASES = "ACGT"


def read_fasta_sequence(path: Path, max_bases: int | None) -> str:
    chunks = []
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip().upper()
            if not line or line.startswith(">"):
                continue
            if max_bases is not None:
                remaining = max_bases - total
                if remaining <= 0:
                    break
                line = line[:remaining]
            chunks.append(line)
            total += len(line)
    return "".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument(
        "--max-bases",
        type=int,
        default=500_000,
        help="Temporary baseline cap per genome. Use 0 for full genomes.",
    )
    args = parser.parse_args()
    max_bases = args.max_bases or None

    metadata = pd.read_csv(METADATA_PATH)
    metadata["genome_id"] = metadata["genome_id"].astype(str)
    kmers = ["".join(parts) for parts in product(BASES, repeat=args.k)]

    genome_ids = []
    corpus = []
    for genome_id in metadata["genome_id"]:
        fasta_path = FASTA_DIR / f"{genome_id}.fna"
        if not fasta_path.exists():
            print(f"Skipping missing FASTA: {genome_id}")
            continue
        print(f"Reading FASTA: {genome_id}")
        genome_ids.append(genome_id)
        corpus.append(read_fasta_sequence(fasta_path, max_bases))

    vectorizer = CountVectorizer(
        analyzer="char",
        ngram_range=(args.k, args.k),
        vocabulary=kmers,
        lowercase=False,
    )
    matrix = vectorizer.transform(corpus)
    totals = matrix.sum(axis=1).A1
    totals[totals == 0] = 1
    normalized = matrix.multiply(1 / totals[:, None]).toarray()

    frame = pd.DataFrame(normalized, columns=[f"kmer_{kmer}" for kmer in kmers])
    frame.insert(0, "genome_id", genome_ids)

    FEATURES_PATH.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(FEATURES_PATH, index=False)
    print(f"Wrote {FEATURES_PATH}")
    print(f"Rows: {len(frame)}")
    print(f"Features: {len(kmers)}")


if __name__ == "__main__":
    main()

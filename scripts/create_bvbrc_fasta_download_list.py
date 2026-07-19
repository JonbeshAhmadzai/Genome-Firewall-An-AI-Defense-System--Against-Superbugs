from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
DEFAULT_GENOME_LIST = ROOT / "data" / "interim" / "bvbrc_genome_ids.txt"
DEFAULT_URL_LIST = ROOT / "data" / "interim" / "bvbrc_fasta_urls.txt"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--genome-list-out", type=Path, default=DEFAULT_GENOME_LIST)
    parser.add_argument("--url-list-out", type=Path, default=DEFAULT_URL_LIST)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    if "genome_id" not in metadata.columns:
        raise SystemExit(f"{args.metadata} does not contain a genome_id column")

    genome_ids = sorted(set(metadata["genome_id"].astype(str)))
    args.genome_list_out.parent.mkdir(parents=True, exist_ok=True)
    args.genome_list_out.write_text("\n".join(genome_ids) + "\n", encoding="utf-8")

    urls = [f"ftps://ftp.bv-brc.org/genomes/{genome_id}/{genome_id}.fna" for genome_id in genome_ids]
    args.url_list_out.write_text("\n".join(urls) + "\n", encoding="utf-8")

    print(f"Genome IDs: {len(genome_ids)}")
    print(f"Wrote {args.genome_list_out}")
    print(f"Wrote {args.url_list_out}")


if __name__ == "__main__":
    main()


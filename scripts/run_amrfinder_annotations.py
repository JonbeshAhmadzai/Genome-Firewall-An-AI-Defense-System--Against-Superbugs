from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FASTA_DIR = ROOT / "data" / "raw" / "fasta"
AMRFINDER_DIR = ROOT / "data" / "interim" / "amrfinder"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", default="Escherichia")
    args = parser.parse_args()

    amrfinder = shutil.which("amrfinder")
    if not amrfinder:
        raise SystemExit(
            "amrfinder is not installed. Install ncbi-amrfinderplus or run this step in WSL/Docker."
        )

    AMRFINDER_DIR.mkdir(parents=True, exist_ok=True)
    fasta_files = sorted(FASTA_DIR.glob("*.fna"))
    if not fasta_files:
        raise SystemExit(f"No FASTA files found in {FASTA_DIR}")

    for fasta_path in fasta_files:
        output_path = AMRFINDER_DIR / f"{fasta_path.stem}.tsv"
        if output_path.exists() and output_path.stat().st_size > 0:
            continue
        command = [
            amrfinder,
            "-n",
            str(fasta_path),
            "-O",
            args.organism,
            "--plus",
            "-o",
            str(output_path),
        ]
        print(" ".join(command))
        subprocess.run(command, check=True)

    print(f"Wrote AMRFinderPlus TSVs to {AMRFINDER_DIR}")


if __name__ == "__main__":
    main()


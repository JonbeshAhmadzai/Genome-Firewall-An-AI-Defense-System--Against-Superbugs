"""Safe, repeatable AMRFinderPlus execution for one assembled genome."""

from __future__ import annotations

import subprocess
import gzip
import shutil
import tempfile
from pathlib import Path


def run_amrfinder(
    genome_fasta: Path,
    output_tsv: Path,
    *,
    executable: str = "amrfinder",
    database: str | None = None,
    organism: str | None = "Escherichia",
    plus: bool = True,
    threads: int = 8,
) -> Path:
    """Annotate ``genome_fasta`` and atomically write the tabular result."""

    genome_fasta = Path(genome_fasta)
    output_tsv = Path(output_tsv)
    if threads < 1:
        raise ValueError("threads must be positive")
    if not genome_fasta.is_file() or genome_fasta.stat().st_size == 0:
        raise FileNotFoundError(f"Genome FASTA is missing or empty: {genome_fasta}")
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    partial = output_tsv.with_suffix(output_tsv.suffix + ".part")
    partial.unlink(missing_ok=True)
    temporary_fasta: Path | None = None
    input_fasta = genome_fasta
    if genome_fasta.suffix == ".gz":
        handle = tempfile.NamedTemporaryFile(prefix="genome-firewall-", suffix=".fna", delete=False)
        temporary_fasta = Path(handle.name)
        handle.close()
        with gzip.open(genome_fasta, "rb") as source, temporary_fasta.open("wb") as destination:
            shutil.copyfileobj(source, destination)
        input_fasta = temporary_fasta
    command = [executable, "-n", str(input_fasta), "-o", str(partial), "--threads", str(threads)]
    if organism:
        command.extend(["--organism", organism])
    if plus:
        command.append("--plus")
    if database:
        command.extend(["--database", database])
    try:
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except FileNotFoundError as error:
            raise RuntimeError(f"AMRFinderPlus executable not found: {executable}") from error
        except subprocess.CalledProcessError as error:
            details = (error.stderr or error.stdout or "").strip()
            raise RuntimeError(f"AMRFinderPlus failed for {genome_fasta}: {details}") from error
    finally:
        if temporary_fasta is not None:
            temporary_fasta.unlink(missing_ok=True)
    if not partial.exists() or partial.stat().st_size == 0:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"AMRFinderPlus produced no output for {genome_fasta}")
    partial.replace(output_tsv)
    return output_tsv

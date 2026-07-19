from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    parts = resolved.parts[1:]
    return "/mnt/" + drive + "/" + "/".join(parts).replace("\\", "/")


def main() -> None:
    project = to_wsl_path(ROOT)
    fasta = f"{project}/data/raw/fasta"
    out = f"{project}/data/interim/amrfinder"

    print("Use Windows conda for Python/modeling. Use WSL only for AMRFinderPlus.")
    print()
    print("1. In an Ubuntu/Debian WSL shell, install micromamba if needed:")
    print(
        "curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | "
        "tar -xvj bin/micromamba"
    )
    print('export MAMBA_ROOT_PREFIX="$HOME/micromamba"')
    print('eval "$($HOME/bin/micromamba shell hook -s bash)"')
    print()
    print("2. Create a tiny WSL env with only AMRFinderPlus:")
    print("micromamba create -n amr -c conda-forge -c bioconda ncbi-amrfinderplus -y")
    print("micromamba activate amr")
    print("amrfinder_update")
    print()
    print("3. Run AMRFinderPlus on the Windows project FASTA files:")
    print(f"cd {project}")
    print(f"mkdir -p {out}")
    print(f"for f in {fasta}/*.fna; do")
    print('  id="$(basename "$f" .fna)"')
    print(f'  if [ ! -s "{out}/${{id}}.tsv" ]; then')
    print(f'    amrfinder -n "$f" -O Escherichia --plus -o "{out}/${{id}}.tsv"')
    print("  fi")
    print("done")
    print()
    print("4. Back in Windows PowerShell, build AMRFinderPlus features and train:")
    print("conda run -n genome python scripts/build_amrfinder_features.py")
    print(
        "conda run -n genome python scripts/train_models.py "
        "--features data/interim/ecoli_amrfinder_features.csv "
        "--feature-prefix amr_ --feature-type amrfinderplus_gene_mutation"
    )


if __name__ == "__main__":
    main()


from __future__ import annotations

import subprocess
import sys
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print(" ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--likely-to-work-max", type=float, default=0.30)
    parser.add_argument("--likely-to-fail-min", type=float, default=0.70)
    args = parser.parse_args()

    if args.likely_to_work_max >= args.likely_to_fail_min:
        raise SystemExit("--likely-to-work-max must be lower than --likely-to-fail-min")

    run([sys.executable, "scripts/build_amrfinder_features.py"])
    run(
        [
            sys.executable,
            "scripts/train_models.py",
            "--features",
            "data/interim/ecoli_amrfinder_features.csv",
            "--feature-prefix",
            "amr_",
            "--feature-type",
            "amrfinderplus_gene_mutation",
            "--output-dir",
            "models/amrfinder",
            "--metrics-out",
            "reports/metrics/amrfinder_model_metrics.csv",
            "--predictions-out",
            "reports/metrics/amrfinder_heldout_predictions.csv",
            "--manifest-out",
            "models/amrfinder/manifest.json",
            "--likely-to-work-max",
            str(args.likely_to_work_max),
            "--likely-to-fail-min",
            str(args.likely_to_fail_min),
        ]
    )


if __name__ == "__main__":
    main()

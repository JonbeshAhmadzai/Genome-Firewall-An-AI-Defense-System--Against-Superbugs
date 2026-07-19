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

    run(
        [
            sys.executable,
            "scripts/train_amrfinder_models.py",
            "--likely-to-work-max",
            str(args.likely_to_work_max),
            "--likely-to-fail-min",
            str(args.likely_to_fail_min),
        ]
    )
    run([sys.executable, "scripts/predict_amrfinder_models.py"])
    run(
        [
            sys.executable,
            "scripts/generate_evaluation_report.py",
            "--metrics",
            "reports/metrics/amrfinder_model_metrics.csv",
            "--predictions",
            "reports/metrics/amrfinder_heldout_predictions.csv",
            "--model-card",
            "reports/amrfinder_model_card.md",
            "--prefix",
            "amrfinder",
            "--feature-source",
            "AMRFinderPlus AMR gene and mutation presence features",
        ]
    )


if __name__ == "__main__":
    main()

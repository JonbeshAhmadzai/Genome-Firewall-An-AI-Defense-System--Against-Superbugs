#!/usr/bin/env python3
"""Validate the config-driven pathogen and target registry."""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from targets_config import ACTIVE_TARGETS, validate_config


def main() -> int:
    errors = validate_config()
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"Configuration is valid. Active bacterial targets: {len(ACTIVE_TARGETS)}")
    for species, antibiotic in ACTIVE_TARGETS:
        print(f"- {species} / {antibiotic}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

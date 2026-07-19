"""Leakage-safe train/calibration/test splitting by genetic group."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit


@dataclass(frozen=True)
class SplitAssignment:
    frame: pd.DataFrame
    group_column: str
    random_state: int

    @property
    def train(self) -> pd.DataFrame:
        return self.frame[self.frame["split"].eq("train")]

    @property
    def calibration(self) -> pd.DataFrame:
        return self.frame[self.frame["split"].eq("calibration")]

    @property
    def test(self) -> pd.DataFrame:
        return self.frame[self.frame["split"].eq("test")]


def _validate_groups(frame: pd.DataFrame, group_column: str) -> None:
    if group_column not in frame.columns:
        raise ValueError(f"Missing required genetic group column {group_column!r}")
    if frame[group_column].isna().any() or frame[group_column].astype(str).str.strip().eq("").any():
        raise ValueError("Every genome must have a non-empty genetic group ID")
    if frame[group_column].nunique() < 3:
        raise ValueError("At least three genetic groups are required for train/calibration/test")


def grouped_three_way_split(
    frame: pd.DataFrame,
    *,
    group_column: str = "homology_group_id",
    test_size: float = 0.20,
    calibration_size: float = 0.20,
    random_state: int = 42,
) -> SplitAssignment:
    """Assign every row to a group-disjoint train/calibration/test split.

    The test split is selected first. Calibration is selected from the
    remaining groups, leaving the rest for training. The returned frame is a
    copy and preserves the caller's row index in ``_original_index``.
    """

    if not 0 < test_size < 1 or not 0 < calibration_size < 1 or test_size + calibration_size >= 1:
        raise ValueError("test_size and calibration_size must be in (0, 1) and sum to less than 1")
    if frame.empty:
        raise ValueError("Cannot split an empty frame")
    _validate_groups(frame, group_column)

    working = frame.copy()
    working["_original_index"] = working.index
    groups = working[group_column].astype(str).to_numpy()
    all_indices = working.index.to_numpy()

    test_splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_cal_idx, test_idx = next(test_splitter.split(all_indices, groups=groups))
    train_cal = working.iloc[train_cal_idx]
    test = working.iloc[test_idx]

    calibration_splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=calibration_size / (1.0 - test_size),
        random_state=random_state,
    )
    train_idx, calibration_idx = next(
        calibration_splitter.split(train_cal.index.to_numpy(), groups=train_cal[group_column].astype(str).to_numpy())
    )
    train = train_cal.iloc[train_idx]
    calibration = train_cal.iloc[calibration_idx]

    result = working.copy()
    result["split"] = ""
    result.loc[train.index, "split"] = "train"
    result.loc[calibration.index, "split"] = "calibration"
    result.loc[test.index, "split"] = "test"
    if result["split"].eq("").any():
        raise RuntimeError("Internal split error: at least one row has no assignment")

    groups_by_split = {
        split: set(result.loc[result["split"].eq(split), group_column].astype(str))
        for split in ("train", "calibration", "test")
    }
    for left in groups_by_split:
        for right in groups_by_split:
            if left < right and groups_by_split[left] & groups_by_split[right]:
                raise RuntimeError(f"Genetic groups overlap between {left} and {right}")

    return SplitAssignment(result, group_column, random_state)


def split_receipt(assignment: SplitAssignment) -> dict[str, object]:
    """Return a JSON-serializable receipt for QA and human leakage review."""

    frame = assignment.frame
    return {
        "group_column": assignment.group_column,
        "random_state": assignment.random_state,
        "rows": {split: int(frame["split"].eq(split).sum()) for split in ("train", "calibration", "test")},
        "groups": {
            split: sorted(frame.loc[frame["split"].eq(split), assignment.group_column].astype(str).unique().tolist())
            for split in ("train", "calibration", "test")
        },
    }

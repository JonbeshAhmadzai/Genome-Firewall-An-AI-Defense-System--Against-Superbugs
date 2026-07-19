from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
AMRFINDER_DIR = ROOT / "data" / "interim" / "amrfinder"
FEATURES_PATH = ROOT / "data" / "interim" / "ecoli_amrfinder_features.csv"
EVIDENCE_PATH = ROOT / "data" / "interim" / "ecoli_amrfinder_evidence.csv"
ALL_HITS_PATH = ROOT / "data" / "interim" / "ecoli_amrfinder_all_hits.csv"


def clean_feature_name(value: object) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "unknown"


def find_first_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    normalized = {column.lower().strip(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]
    return None


def feature_from_row(row: pd.Series, symbol_column: str, subtype_column: str | None) -> str:
    symbol = clean_feature_name(row.get(symbol_column))
    subtype = clean_feature_name(row.get(subtype_column)) if subtype_column else ""
    if subtype and subtype != symbol:
        return f"amr_{symbol}_{subtype}"
    return f"amr_{symbol}"


def get_value(row: pd.Series, column: str | None) -> object:
    if not column:
        return ""
    return row.get(column, "")


def main() -> None:
    tsv_files = sorted(AMRFINDER_DIR.glob("*.tsv"))
    if not tsv_files:
        raise SystemExit(
            f"No AMRFinderPlus TSV files found in {AMRFINDER_DIR}. Run scripts/run_amrfinder_annotations.py first."
        )

    rows = []
    evidence_rows = []
    all_hit_rows = []
    all_features = set()
    for tsv_path in tsv_files:
        genome_id = tsv_path.stem
        frame = pd.read_csv(tsv_path, sep="\t")
        if frame.empty:
            rows.append({"genome_id": genome_id})
            continue

        symbol_column = find_first_column(frame, ["Element symbol", "Gene symbol", "Gene", "Element"])
        subtype_column = find_first_column(frame, ["Element subtype", "Subtype", "Mutation"])
        type_column = find_first_column(frame, ["Type"])
        class_column = find_first_column(frame, ["Class", "Subclass"])
        subclass_column = find_first_column(frame, ["Subclass"])
        method_column = find_first_column(frame, ["Method", "Element type"])
        name_column = find_first_column(frame, ["Element name", "Closest reference name"])
        identity_column = find_first_column(frame, ["% Identity to reference"])
        coverage_column = find_first_column(frame, ["% Coverage of reference"])
        if not symbol_column:
            raise ValueError(f"Could not identify AMRFinderPlus element column in {tsv_path}")

        feature_values = {}
        for _, row in frame.iterrows():
            hit_type = str(get_value(row, type_column)).strip().upper()
            hit_class = str(get_value(row, class_column)).strip()
            hit_subclass = str(get_value(row, subclass_column)).strip()
            hit = {
                "genome_id": genome_id,
                "element_symbol": get_value(row, symbol_column),
                "element_name": get_value(row, name_column),
                "type": hit_type,
                "subtype": get_value(row, subtype_column),
                "class": hit_class,
                "subclass": hit_subclass,
                "method": get_value(row, method_column),
                "identity_to_reference": get_value(row, identity_column),
                "coverage_of_reference": get_value(row, coverage_column),
            }
            all_hit_rows.append(hit)

            if hit_type != "AMR":
                continue

            feature = feature_from_row(row, symbol_column, subtype_column)
            feature_values[feature] = 1
            all_features.add(feature)
            evidence_rows.append(
                {
                    "genome_id": genome_id,
                    "feature": feature,
                    **hit,
                }
            )

        rows.append({"genome_id": genome_id, **feature_values})

    features = pd.DataFrame(rows).fillna(0)
    for feature in sorted(all_features):
        if feature not in features.columns:
            features[feature] = 0
    ordered_columns = ["genome_id"] + sorted(column for column in features.columns if column != "genome_id")
    features = features[ordered_columns]

    FEATURES_PATH.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(FEATURES_PATH, index=False)
    pd.DataFrame(evidence_rows).to_csv(EVIDENCE_PATH, index=False)
    pd.DataFrame(all_hit_rows).to_csv(ALL_HITS_PATH, index=False)
    print(f"Wrote {FEATURES_PATH}")
    print(f"Wrote {EVIDENCE_PATH}")
    print(f"Wrote {ALL_HITS_PATH}")
    print(f"Genomes: {len(features)}")
    print(f"AMR features: {len(ordered_columns) - 1}")


if __name__ == "__main__":
    main()

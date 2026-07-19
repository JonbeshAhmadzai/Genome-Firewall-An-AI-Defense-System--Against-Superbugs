from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

METADATA_PATH = ROOT / "data" / "processed" / "ecoli_genome_metadata.csv"
FEATURES_PATH = ROOT / "data" / "interim" / "bvbrc_genome_features.csv"
ANTIBIOTICS_PATH = ROOT / "configs" / "antibiotics.yaml"
OUT_PATH = ROOT / "data" / "interim" / "target_presence.csv"


def normalize_text(value: object) -> str:
    text = str(value or "").lower()
    return re.sub(r"\s+", " ", text).strip()


def normalize_gene(value: object) -> str:
    text = normalize_text(value)
    return re.sub(r"[^a-z0-9_/-]+", "", text)


def find_target_matches(features: pd.DataFrame, marker_config: dict) -> list[str]:
    genes = {normalize_gene(value) for value in marker_config.get("genes", [])}
    keywords = [normalize_text(value) for value in marker_config.get("product_keywords", [])]
    matches: set[str] = set()

    if genes and "_gene_norm" in features.columns:
        observed_genes = set(features["_gene_norm"].dropna())
        for gene in sorted(genes & observed_genes):
            matches.add(f"gene:{gene}")
    if keywords and "_product_norm" in features.columns:
        product_text = features["_product_norm"].dropna()
        for keyword in keywords:
            if keyword and product_text.str.contains(re.escape(keyword), regex=True, na=False).any():
                matches.add(f"product:{keyword}")

    return sorted(matches)


def build_target_presence(metadata: pd.DataFrame, features: pd.DataFrame, antibiotics: dict) -> pd.DataFrame:
    metadata = metadata.copy()
    metadata["genome_id"] = metadata["genome_id"].astype(str)
    features = features.copy()
    if features.empty:
        features = pd.DataFrame(columns=["genome_id"])
    features["genome_id"] = features["genome_id"].astype(str)
    if "fetch_error" not in features.columns:
        features["fetch_error"] = ""
    if "gene" not in features.columns:
        features["gene"] = ""
    if "product" not in features.columns:
        features["product"] = ""
    features["_gene_norm"] = features["gene"].map(normalize_gene)
    features["_product_norm"] = features["product"].map(normalize_text)

    feature_groups = {genome_id: group for genome_id, group in features.groupby("genome_id")}
    rows: list[dict] = []
    for genome_id in metadata["genome_id"].dropna().astype(str).unique():
        genome_features = feature_groups.get(genome_id, pd.DataFrame())
        valid_features = genome_features[
            genome_features.get("fetch_error", pd.Series("", index=genome_features.index)).fillna("").eq("")
        ].copy()
        fetch_errors = sorted(
            {
                str(value)
                for value in genome_features.get("fetch_error", pd.Series(dtype=str)).dropna()
                if str(value).strip()
            }
        )

        for antibiotic, drug_config in antibiotics.get("antibiotics", {}).items():
            markers = drug_config.get("target_markers", {})
            targets = ", ".join(drug_config.get("primary_targets", []))
            sources = drug_config.get("target_sources", [])
            if valid_features.empty:
                status = "unknown_annotation"
                matched = []
                reason = "; ".join(fetch_errors) if fetch_errors else "no genome features available"
            else:
                matched = find_target_matches(valid_features, markers)
                status = "present" if matched else "absent"
                reason = "matched configured target marker" if matched else "no configured target marker matched"

            rows.append(
                {
                    "genome_id": genome_id,
                    "antibiotic": antibiotic,
                    "target_status": status,
                    "target_present": status == "present",
                    "matched_target_markers": "; ".join(matched[:12]),
                    "molecular_targets": targets,
                    "target_confidence": drug_config.get("target_confidence", "not_configured"),
                    "target_source_count": len(sources),
                    "target_gate_reason": reason,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build per-genome per-antibiotic target presence from BV-BRC genome features."
    )
    parser.add_argument("--metadata", type=Path, default=METADATA_PATH)
    parser.add_argument("--features", type=Path, default=FEATURES_PATH)
    parser.add_argument("--antibiotics", type=Path, default=ANTIBIOTICS_PATH)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    if args.features.exists():
        features = pd.read_csv(args.features, dtype=str)
    else:
        features = pd.DataFrame(columns=["genome_id", "fetch_error"])
    antibiotics = yaml.safe_load(args.antibiotics.read_text(encoding="utf-8"))
    target_presence = build_target_presence(metadata, features, antibiotics)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    target_presence.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(target_presence)} rows", flush=True)
    if not target_presence.empty:
        print(
            target_presence.pivot_table(
                index="antibiotic",
                columns="target_status",
                values="genome_id",
                aggfunc="count",
                fill_value=0,
            ).to_string(),
            flush=True,
        )


if __name__ == "__main__":
    main()

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "dataset.yaml"
RAW_PATH = ROOT / "data" / "raw" / "bvbrc_ecoli_lab_amr.jsonl"
LABELS_PATH = ROOT / "data" / "processed" / "bvbrc_ecoli_labels.csv"
COUNTS_PATH = ROOT / "reports" / "metrics" / "bvbrc_ecoli_label_counts.csv"


VALID_NAME_RE = re.compile(r"^[a-z][a-z0-9 +./()-]{1,80}$")


def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def normalize_label(value: object, label_values: dict) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip().lower()
    for label, accepted in label_values.items():
        if normalized in accepted:
            return label
    return None


def valid_antibiotic(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip().lower()
    if not VALID_NAME_RE.match(text):
        return False
    blocked = {"publication", "publicatiion", "unknown", "not defined", "not specified"}
    return text not in blocked


def fetch_rows(config: dict) -> list[dict]:
    base = config["api_base_url"].rstrip("/")
    page_size = int(config["page_size"])
    max_rows = int(config["max_rows"])
    species = quote(config["species"])
    evidence = quote(config["evidence"])

    query = f"and(keyword({species}),eq(evidence,{evidence}))"
    fields = ",".join(
        [
            "genome_id",
            "genome_name",
            "taxon_id",
            "antibiotic",
            "resistant_phenotype",
            "measurement",
            "measurement_value",
            "measurement_unit",
            "measurement_sign",
            "evidence",
            "laboratory_typing_method",
            "laboratory_typing_platform",
        ]
    )

    rows: list[dict] = []
    session = requests.Session()
    for start in range(0, max_rows, page_size):
        url = (
            f"{base}/genome_amr/?{query}"
            f"&select({fields})&limit({page_size},{start})"
            "&http_accept=application/json"
        )
        response = session.get(url, timeout=60)
        response.raise_for_status()
        batch = response.json()
        if not batch:
            break
        rows.extend(batch)
        print(f"Fetched {len(rows)} rows")
        if len(batch) < page_size:
            break
    return rows


def clean_rows(rows: list[dict], config: dict) -> pd.DataFrame:
    clean: list[dict] = []
    for row in rows:
        antibiotic = str(row.get("antibiotic") or "").strip().lower()
        label = normalize_label(row.get("resistant_phenotype"), config["label_values"])
        if not label or label == "uncertain":
            continue
        if not valid_antibiotic(antibiotic):
            continue
        clean.append(
            {
                "genome_id": row.get("genome_id"),
                "genome_name": row.get("genome_name"),
                "taxon_id": row.get("taxon_id"),
                "antibiotic": antibiotic,
                "label": label,
                "resistant_phenotype": row.get("resistant_phenotype"),
                "measurement": row.get("measurement"),
                "measurement_value": row.get("measurement_value"),
                "measurement_unit": row.get("measurement_unit"),
                "measurement_sign": row.get("measurement_sign"),
                "evidence": row.get("evidence"),
                "laboratory_typing_method": row.get("laboratory_typing_method"),
                "laboratory_typing_platform": row.get("laboratory_typing_platform"),
            }
        )

    frame = pd.DataFrame(clean)
    if frame.empty:
        return frame

    frame = frame.drop_duplicates(["genome_id", "antibiotic", "label"])
    return frame.sort_values(["antibiotic", "genome_id"]).reset_index(drop=True)


def write_outputs(raw_rows: list[dict], labels: pd.DataFrame) -> None:
    RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    LABELS_PATH.parent.mkdir(parents=True, exist_ok=True)
    COUNTS_PATH.parent.mkdir(parents=True, exist_ok=True)

    with RAW_PATH.open("w", encoding="utf-8") as handle:
        for row in raw_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    labels.to_csv(LABELS_PATH, index=False)

    if labels.empty:
        COUNTS_PATH.write_text("antibiotic,label,count\n", encoding="utf-8")
        return

    counts = (
        labels.groupby(["antibiotic", "label"])
        .size()
        .reset_index(name="count")
        .sort_values(["antibiotic", "label"])
    )
    counts.to_csv(COUNTS_PATH, index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> None:
    config = load_config()
    raw_rows = fetch_rows(config)
    labels = clean_rows(raw_rows, config)
    write_outputs(raw_rows, labels)

    print(f"Raw rows: {len(raw_rows)}")
    print(f"Usable labeled rows: {len(labels)}")
    if not labels.empty:
        print(f"Genomes: {labels['genome_id'].nunique()}")
        print(f"Antibiotics: {labels['antibiotic'].nunique()}")
        print("Top antibiotics by usable labels:")
        print(Counter(labels["antibiotic"]).most_common(10))
    print(f"Wrote {LABELS_PATH}")


if __name__ == "__main__":
    main()

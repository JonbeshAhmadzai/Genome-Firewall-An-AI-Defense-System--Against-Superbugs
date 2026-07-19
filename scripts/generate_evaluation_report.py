from __future__ import annotations

import sys
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from sklearn.calibration import calibration_curve
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

FIGURES_DIR = ROOT / "reports" / "figures"
ANTIBIOTICS_PATH = ROOT / "configs" / "antibiotics.yaml"


def decision_to_label(decision: str) -> int | None:
    if decision == "likely_to_fail":
        return 1
    if decision == "likely_to_work":
        return 0
    return None


def write_reliability_plots(predictions: pd.DataFrame, prefix: str) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for antibiotic, group in predictions.groupby("antibiotic"):
        y_true = group["true_label"].map({"susceptible": 0, "resistant": 1}).to_numpy()
        y_prob = group["resistant_probability"].to_numpy()
        if len(set(y_true)) < 2:
            continue
        prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=min(5, len(group)), strategy="uniform")
        plt.figure(figsize=(5, 4))
        plt.plot([0, 1], [0, 1], linestyle="--", color="0.5", label="ideal")
        plt.plot(prob_pred, prob_true, marker="o", label=antibiotic)
        plt.xlabel("Mean predicted resistant probability")
        plt.ylabel("Observed resistant fraction")
        plt.title(f"Reliability: {antibiotic}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / f"{prefix}_reliability_{antibiotic}.png", dpi=160)
        plt.close()


def write_decision_summaries(predictions: pd.DataFrame, prefix: str) -> None:
    summary = (
        predictions.groupby(["antibiotic", "decision"])
        .size()
        .reset_index(name="count")
        .sort_values(["antibiotic", "decision"])
    )
    summary.to_csv(ROOT / "reports" / "metrics" / f"{prefix}_decision_summary.csv", index=False)

    for antibiotic, group in predictions.groupby("antibiotic"):
        called = group[group["decision"] != "no_call"].copy()
        if called.empty:
            continue
        called["y_true"] = called["true_label"].map({"susceptible": 0, "resistant": 1})
        called["y_pred"] = called["decision"].map(decision_to_label)
        matrix = confusion_matrix(called["y_true"], called["y_pred"], labels=[0, 1])
        display = ConfusionMatrixDisplay(matrix, display_labels=["work", "fail"])
        display.plot(values_format="d", colorbar=False)
        plt.title(f"Called Decision Matrix: {antibiotic}")
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / f"{prefix}_called_confusion_{antibiotic}.png", dpi=160)
        plt.close()


def write_model_card(metrics: pd.DataFrame, predictions: pd.DataFrame, model_card_path: Path, feature_source: str) -> None:
    antibiotic_config = yaml.safe_load(ANTIBIOTICS_PATH.read_text(encoding="utf-8"))
    supported = ", ".join(metrics["antibiotic"].tolist())
    no_call = (
        predictions.groupby("antibiotic")["decision"]
        .apply(lambda values: float((values == "no_call").mean()))
        .reset_index(name="no_call_rate")
    )

    metric_columns = [
        "antibiotic",
        "balanced_accuracy",
        "resistant_recall",
        "susceptible_recall",
        "f1",
        "auroc",
        "pr_auc",
        "brier_score",
    ]
    no_call_columns = ["antibiotic", "no_call_rate"]

    if "amrfinder" in feature_source.lower():
        feature_limitation = "- AMRFinderPlus annotations are used as features; regenerate this card after all TSVs finish so metrics reflect the full synchronized cohort."
        explanation_line = "- Honest explanations: drug-specific AMRFinderPlus hits are reported separately from statistical model associations; statistical associations are not proof of biological causality."
    else:
        feature_limitation = "- Current scores are not submission-quality because the temporary baseline is not AMRFinderPlus-based."
        explanation_line = "- Honest explanations: current k-mer evidence is marked as statistical association only; AMRFinderPlus features will support known gene/mutation evidence."

    lines = [
        "# Genome Firewall Model Card",
        "",
        "## Intended Use",
        "",
        "Research prototype for defensive antibiotic-response decision support from reconstructed bacterial genomes. All outputs require confirmation by standard laboratory testing.",
        "",
        "## Current Scope",
        "",
        f"- Species: `{antibiotic_config['supported_species']}`",
        f"- Trained antibiotics: {supported}",
        f"- Current feature source: {feature_source}.",
        "- Unsupported: organism design, modification, optimization, treatment automation, raw sample processing, species identification.",
        "",
        "## Decision Policy",
        "",
        "- `likely_to_fail`: resistant probability >= 0.70",
        "- `likely_to_work`: resistant probability <= 0.30",
        "- `no_call`: probability between those thresholds or insufficient evidence.",
        "",
        "## Held-Out Metrics",
        "",
        markdown_table(metrics[metric_columns]),
        "",
        "## No-Call Rates",
        "",
        markdown_table(no_call[no_call_columns]),
        "",
        "## Responsibility Requirements",
        "",
        "- Defensive by construction: predicts resistance that may already exist; does not generate or suggest organism changes.",
        "- Honest generalization: grouped split uses `cgmlst_hc100` when available; training writes a split-audit CSV with train/test group overlap counts.",
        "- Calibrated confidence: model uses sigmoid calibration and reports Brier score plus reliability plots.",
        "- No-call option: uncertainty is surfaced as a first-class output.",
        explanation_line,
        "- Human oversight: every report includes standard-lab-confirmation language.",
        "",
        "## Limitations",
        "",
        feature_limitation,
        "- Confidence and generalization estimates are unstable until the AMRFinderPlus run covers the synchronized FASTA-backed cohort.",
        "- Target gate verifies configured target markers from BV-BRC genome annotations when `data/interim/target_presence.csv` is refreshed for the current cohort; uploaded FASTAs still need an annotation path.",
    ]
    model_card_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def markdown_table(frame: pd.DataFrame) -> str:
    display = frame.copy()
    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    header = "| " + " | ".join(display.columns) + " |"
    separator = "| " + " | ".join("---" for _ in display.columns) + " |"
    rows = ["| " + " | ".join(str(value) for value in row) + " |" for row in display.to_numpy()]
    return "\n".join([header, separator, *rows])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", type=Path, default=ROOT / "reports" / "metrics" / "model_metrics.csv")
    parser.add_argument("--predictions", type=Path, default=ROOT / "reports" / "metrics" / "heldout_predictions.csv")
    parser.add_argument("--model-card", type=Path, default=ROOT / "reports" / "model_card.md")
    parser.add_argument("--prefix", default="kmer")
    parser.add_argument("--feature-source", default="temporary capped 4-mer frequencies")
    args = parser.parse_args()

    metrics = pd.read_csv(args.metrics)
    predictions = pd.read_csv(args.predictions)
    write_reliability_plots(predictions, args.prefix)
    write_decision_summaries(predictions, args.prefix)
    write_model_card(metrics, predictions, args.model_card, args.feature_source)
    print(f"Wrote {args.model_card}")
    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()

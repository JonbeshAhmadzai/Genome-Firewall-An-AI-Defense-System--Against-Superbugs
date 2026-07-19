"""Safe optional LLM explanations for Genome Firewall outputs.

The language model is deliberately kept downstream of the scientific pipeline:
it receives structured prediction/QC/evidence summaries, but it cannot change
the prediction, confidence, target gate, or no-call decision.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from dotenv import load_dotenv
except ImportError:  # pragma: no cover - requirements install provides python-dotenv
    def load_dotenv(*_: Any, **__: Any) -> bool:
        """Fallback that keeps the deterministic pipeline usable without extras."""

        return False


DEFAULT_MODEL = "gpt-4o-mini"
MAX_EVIDENCE_ROWS = 40
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _load_environment() -> None:
    """Load the local .env without overwriting explicitly exported variables."""

    load_dotenv(PROJECT_ROOT / ".env", override=False)


def _configured_api_key() -> str | None:
    key = os.getenv("OPENAI_API_KEY", "").strip()
    if not key or key in {"your-key-here", "replace-with-your-key"}:
        return None
    return key


def llm_available() -> bool:
    """Return whether an API key is configured and LLM explanations are enabled."""

    _load_environment()
    return _configured_api_key() is not None and os.getenv("OPENAI_LLM_ENABLED", "true").lower() not in {
        "0",
        "false",
        "no",
        "off",
    }


def _json_safe(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def build_explanation_payload(
    predictions: pd.DataFrame,
    *,
    species: str,
    antibiotic: str,
    qc: dict[str, Any] | None = None,
    evidence: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Build a bounded, structured payload; never include FASTA or sequences."""

    prediction_columns = [
        column
        for column in (
            "genome_id",
            "antibiotic",
            "prediction",
            "prediction_before_target_gate",
            "probability_resistant",
            "confidence",
            "target_status",
            "target_gate_reason",
            "evidence_category",
        )
        if column in predictions.columns
    ]
    prediction_rows = [
        {column: _json_safe(row[column]) for column in prediction_columns}
        for _, row in predictions.head(100).iterrows()
    ]

    evidence_rows: list[dict[str, Any]] = []
    if evidence is not None and not evidence.empty:
        evidence_columns = [
            column
            for column in (
                "genome_id",
                "gene_symbol",
                "element_name",
                "amr_class",
                "amr_subclass",
                "subtype",
                "evidence",
            )
            if column in evidence.columns
        ]
        evidence_rows = [
            {column: _json_safe(row[column]) for column in evidence_columns}
            for _, row in evidence.head(MAX_EVIDENCE_ROWS).iterrows()
        ]

    return {
        "species": species,
        "antibiotic": antibiotic,
        "qc": {
            key: _json_safe(value)
            for key, value in (qc or {}).items()
            if key not in {"path", "input_path", "temporary_path"}
        },
        "predictions": prediction_rows,
        "amr_evidence": evidence_rows,
    }


def explain_predictions(
    predictions: pd.DataFrame,
    *,
    species: str,
    antibiotic: str,
    qc: dict[str, Any] | None = None,
    evidence: pd.DataFrame | None = None,
) -> str:
    """Generate a plain-language explanation without changing model outputs.

    Raises a clear RuntimeError when the optional dependency/key is unavailable.
    The caller should display a deterministic report in that case.
    """

    _load_environment()
    api_key = _configured_api_key()
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not configured. Add it to .env or the shell environment.")
    if os.getenv("OPENAI_LLM_ENABLED", "true").lower() in {"0", "false", "no", "off"}:
        raise RuntimeError("LLM explanations are disabled by OPENAI_LLM_ENABLED.")

    try:
        from openai import OpenAI
    except ImportError as error:  # pragma: no cover - depends on optional environment
        raise RuntimeError("The openai package is not installed. Run: pip install -r requirements.txt") from error

    payload = build_explanation_payload(
        predictions,
        species=species,
        antibiotic=antibiotic,
        qc=qc,
        evidence=evidence,
    )
    instructions = (
        "You are the explanation layer for Genome Firewall, a research prototype. "
        "Do not act as a predictor and do not recommend, select, or change treatment. "
        "Use only the supplied structured data; never invent missing genes, mutations, "
        "phenotypes, or certainty. Explain each prediction, confidence, target status, "
        "and no-call reason in plain language. If prediction is 'no-call' but "
        "prediction_before_target_gate is a different call, explicitly say that the "
        "target-compatibility gate overrode the model call. A high model probability "
        "does not override an uncertain target gate. Clearly separate detected AMR evidence "
        "from statistical model output. If no AMR evidence is supplied, say that no "
        "reportable determinant was detected and that this does not prove susceptibility. "
        "End with: 'Confirm with standard laboratory AST; "
        "this research prototype is not a treatment decision.'"
    )
    response = OpenAI(api_key=api_key).responses.create(
        model=os.getenv("OPENAI_MODEL", DEFAULT_MODEL),
        instructions=instructions,
        input=json.dumps(payload, ensure_ascii=False),
        store=False,
    )
    text = getattr(response, "output_text", None)
    if not text:
        raise RuntimeError("OpenAI returned no explanation text.")
    return text.strip()

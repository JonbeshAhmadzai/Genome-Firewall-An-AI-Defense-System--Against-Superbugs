"""FastAPI backend for the parallel Genome Firewall HTML demo."""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from src.explanation.llm_explainer import explain_predictions, llm_available
from src.genome_reader.build_features import build_feature_matrix
from src.genome_reader.fasta_qc import fasta_stats
from src.genome_reader.run_amrfinder import run_amrfinder
from src.pipeline.router import resolve_route
from src.predictor.predict import load_model, predict_target
from targets_config import enabled_species


ROOT = Path(__file__).resolve().parent
WEB_DIR = ROOT / "web"
MODEL_ROOT = ROOT / "models" / "cohort100_test25"
METRICS_PATH = ROOT / "reports" / "cohort100_test25" / "model_metrics.csv"
BENCHMARK_PATH = ROOT / "reports" / "cohort100_test25" / "model_benchmark.csv"
FEATURES_PATH = ROOT / "data" / "processed" / "cohort100" / "amrfinder_features.csv"
EVIDENCE_PATH = ROOT / "data" / "processed" / "cohort100" / "amrfinder_evidence.csv"
AMRFINDER = os.getenv("AMRFINDER_EXECUTABLE", "/home/becode/miniconda3/envs/genome-firewall/bin/amrfinder")

app = FastAPI(title="Genome Firewall API", version="0.1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)
app.mount("/static", StaticFiles(directory=WEB_DIR), name="static")


def _metrics() -> pd.DataFrame:
    return pd.read_csv(METRICS_PATH) if METRICS_PATH.exists() else pd.DataFrame()


def _csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _available_species() -> list[str]:
    metrics = _metrics()
    if "species" in metrics and not metrics.empty:
        return sorted(metrics["species"].dropna().astype(str).unique())
    return list(enabled_species(kind="bacterium"))


def _available_drugs(species: str) -> list[str]:
    metrics = _metrics()
    if "antibiotic" not in metrics:
        return []
    rows = metrics
    if "species" in metrics:
        rows = metrics[metrics["species"].astype(str).eq(species)]
    return sorted(rows["antibiotic"].dropna().astype(str).unique())


def _json_records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    if frame.empty:
        return []
    return json.loads(frame.to_json(orient="records"))


@app.get("/", response_class=FileResponse)
def index() -> FileResponse:
    return FileResponse(WEB_DIR / "index.html")


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {"status": "ok", "service": "genome-firewall-api", "llm_available": llm_available()}


@app.get("/api/config")
def config() -> dict[str, Any]:
    species = _available_species()
    return {
        "species": species,
        "drugs": {name: _available_drugs(name) for name in species},
        "llm_available": llm_available(),
        "llm_model": os.getenv("OPENAI_MODEL", "gpt-4o-mini"),
    }


@app.get("/api/overview")
def overview() -> dict[str, Any]:
    """Return the non-sensitive cohort, validation, and comparison metadata."""

    metrics = _metrics()
    benchmark = _csv(BENCHMARK_PATH)
    features = _csv(FEATURES_PATH)
    evidence = _csv(EVIDENCE_PATH)
    genome_count = int(features["genome_id"].nunique()) if "genome_id" in features else 0
    return {
        "summary": {
            "genomes": genome_count,
            "drug_models": int(metrics["antibiotic"].nunique()) if "antibiotic" in metrics else 0,
            "amr_features": max(len(features.columns) - 1, 0) if not features.empty else 0,
            "evidence_rows": len(evidence),
        },
        "scope": {
            "species": _available_species(),
            "antibiotics": sorted(metrics["antibiotic"].dropna().astype(str).unique()) if "antibiotic" in metrics else [],
            "status": "100-genome research MVP; not clinically validated",
        },
        "decision_policy": {
            "minimum_confidence": 0.70,
            "no_call": "Low confidence or unverified target compatibility",
            "safety": "Confirm every result with standard laboratory AST.",
        },
        "validation_split": {
            "train": 0.60,
            "calibration": 0.15,
            "test": 0.25,
            "method": "Group-disjoint split by homology group",
        },
        "pipeline": [
            "FASTA quality control",
            "AMRFinderPlus AMR evidence",
            "Engineered determinant features",
            "Grouped calibrated model prediction",
            "Target compatibility gate",
            "Optional evidence-grounded GPT explanation",
        ],
        "held_out_metrics": _json_records(metrics),
        "model_comparison": _json_records(benchmark),
    }


@app.post("/api/predict")
async def predict(
    file: UploadFile = File(...),
    species: str = Form(...),
    explain: bool = Form(False),
) -> dict[str, Any]:
    """Run FASTA QC, AMRFinderPlus, model scoring, and optional explanation."""

    if species not in _available_species():
        raise HTTPException(status_code=400, detail=f"No trained model scope is available for {species}.")
    drugs = _available_drugs(species)
    if not drugs:
        raise HTTPException(status_code=400, detail=f"No trained antibiotic models are available for {species}.")
    try:
        with tempfile.TemporaryDirectory(prefix="genome-firewall-api-") as temporary:
            safe_name = Path(file.filename or "upload.fna").name
            input_path = Path(temporary) / safe_name
            output_path = Path(temporary) / "amrfinder.tsv"
            input_path.write_bytes(await file.read())
            qc = fasta_stats(input_path)
            route = resolve_route(species, drugs[0], model_root=MODEL_ROOT)
            if not route.supported:
                raise HTTPException(status_code=400, detail=route.unsupported_reason())
            run_amrfinder(
                input_path,
                output_path,
                executable=AMRFINDER,
                organism=route.amrfinder_organism,
                plus=True,
                threads=8,
            )
            matrix, evidence = build_feature_matrix([output_path])
            scored_frames: list[pd.DataFrame] = []
            for drug in drugs:
                route = resolve_route(species, drug, model_root=MODEL_ROOT)
                if not route.supported or not route.model_path.exists():
                    continue
                artifact = load_model(route.model_path)
                for column in artifact["feature_columns"]:
                    if column not in matrix.columns:
                        matrix[column] = 0
                scored = predict_target(
                    artifact,
                    matrix,
                    target_present=None,
                    target_names=route.molecular_targets,
                )
                scored.insert(1, "species", species)
                scored.insert(2, "antibiotic", drug)
                scored["evidence_category"] = (
                    "known resistance gene or DNA change detected"
                    if not evidence.empty
                    else "no known resistance signal found"
                )
                scored_frames.append(scored)
            scored_frame = pd.concat(scored_frames, ignore_index=True) if scored_frames else pd.DataFrame()
            response: dict[str, Any] = {
                "qc": qc,
                "species": species,
                "predictions": _json_records(scored_frame),
                "amr_evidence": _json_records(evidence),
                "safety_warning": "Research prototype only. Confirm with standard laboratory testing.",
                "llm_explanation": None,
            }
            if explain:
                if not llm_available():
                    response["llm_explanation"] = "LLM explanation is unavailable; deterministic results are still returned."
                else:
                    response["llm_explanation"] = explain_predictions(
                        scored_frame,
                        species=species,
                        antibiotic="all supported antibiotics",
                        qc=qc,
                        evidence=evidence,
                    )
            return response
    except HTTPException:
        raise
    except Exception as error:
        raise HTTPException(status_code=422, detail=str(error)) from error

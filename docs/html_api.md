# FastAPI + HTML demo

Streamlit remains available. The parallel browser UI runs with:

```bash
uvicorn api:app --reload --port 8001
```

Open `http://127.0.0.1:8001/`.

Endpoints:

- `GET /api/health` — service and optional LLM status
- `GET /api/config` — trained species and antibiotic choices
- `POST /api/predict` — multipart FASTA upload; returns QC, predictions, AMR evidence, safety warning, and optional explanation

Example request:

```bash
curl -F species='Escherichia coli' \
     -F explain=false \
     -F file=@genome.fna.gz \
     http://127.0.0.1:8001/api/predict
```

The endpoint never sends the FASTA itself to the LLM. Only structured QC,
prediction, and AMR-evidence summaries are eligible for explanation.

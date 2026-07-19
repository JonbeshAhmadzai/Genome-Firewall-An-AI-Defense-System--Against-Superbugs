# Render deployment

The repository includes a Docker-based Render service definition. Docker is
used because the upload endpoint runs AMRFinderPlus, which is a Conda/bioconda
binary rather than a PyPI package.

## Deploy from the Render dashboard

1. Push the branch containing `Dockerfile`, `environment.render.yml`,
   `requirements.render.txt`, `.dockerignore`, and `render.yaml` to GitHub.
2. In Render, choose **New → Blueprint** and select the repository/branch.
3. Render reads `render.yaml`, builds the Docker image, and starts Uvicorn on
   Render's `$PORT`.
4. Open the generated URL and verify `/api/health` returns `status: ok`.
5. If explanations are needed, add `OPENAI_API_KEY` in Render's **Environment**
   settings as a secret. Do not commit the key or put it in `render.yaml`.

The default deployment keeps LLM explanations disabled, so the public demo
does not incur API usage until the key is explicitly configured. The model
artifacts under `models/cohort100_test25/` and the matching reports must remain
tracked; raw genomes and local processed data are deliberately excluded from
the image.

## CLI alternative

Create a Render Web Service from the repository, select **Docker** as the
runtime, and use the defaults from `Dockerfile`. Set the health check path to
`/api/health`. Render automatically supplies `PORT`; do not hard-code it.

## Start commands

Local HTML/FastAPI dashboard:

```bash
uvicorn api:app --host 127.0.0.1 --port 8001 --reload
```

Render's Docker start command (already configured in `Dockerfile`):

```bash
uvicorn api:app --host 0.0.0.0 --port ${PORT:-8000}
```

Optional local Streamlit fallback:

```bash
streamlit run app.py --server.address 127.0.0.1 --server.port 8502
```

## Important MVP limitation

This is a research demo, not a clinical service. Keep the safety warning and
AST confirmation message in the UI. Render's free instance may sleep between
requests, and the first AMRFinderPlus request can be slower while the service
wakes up.

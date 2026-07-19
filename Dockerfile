# Render runtime for the FastAPI + HTML demo.
# Conda supplies AMRFinderPlus; pip supplies the Python API/runtime packages.
FROM continuumio/miniconda3:24.9.2-0

WORKDIR /app

# Install the reproducible runtime before copying source files so Docker can
# reuse this expensive layer when only application code changes.
COPY environment.render.yml requirements.render.txt ./
RUN conda env create -f environment.render.yml \
    && /opt/conda/envs/genome-firewall-render/bin/amrfinder -u \
    && conda clean --all --yes \
    && rm -f environment.render.yml requirements.render.txt

ENV PATH="/opt/conda/envs/genome-firewall-render/bin:${PATH}" \
    PYTHONUNBUFFERED="1" \
    PYTHONDONTWRITEBYTECODE="1" \
    AMRFINDER_EXECUTABLE="amrfinder"

# .dockerignore keeps secrets and local/raw data out of the image.
COPY . .

EXPOSE 8000
CMD ["sh", "-c", "uvicorn api:app --host 0.0.0.0 --port ${PORT:-8000}"]

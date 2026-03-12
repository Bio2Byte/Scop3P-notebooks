# syntax=docker/dockerfile:1.7
FROM python:3.12-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ && \
    rm -rf /var/lib/apt/lists/* && \
    pip install -U pip

COPY requirements-biophysics.txt /app/requirements-biophysics.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install -r /app/requirements-biophysics.txt

COPY requirements-shiny.txt /app/requirements-shiny.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install -r /app/requirements-shiny.txt

COPY --chmod=755 TM-align /usr/local/bin/TM-align
COPY . /app

EXPOSE 8000

FROM runtime AS peptide-mapper
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/app/apps/peptide_mapper/app.py"]

FROM runtime AS structure-viz
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/app/apps/structure_viz/app.py"]

FROM runtime AS mutation-effect
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/app/apps/mutation_effect/app.py"]

FROM runtime AS scop3p-suite
CMD ["python", "-m", "uvicorn", "apps.portal.main:app", "--host", "0.0.0.0", "--port", "8000"]

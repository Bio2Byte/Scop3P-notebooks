# syntax=docker/dockerfile:1.7
FROM python:3.12-slim AS runtime

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION=dev

LABEL org.opencontainers.image.title="Scop3P-Toolkit" \
    org.opencontainers.image.description="Scop3P-Toolkit: Tools for exploring and extending Scop3P" \
    org.opencontainers.image.url="https://github.com/Bio2Byte/Scop3P-notebooks" \
    org.opencontainers.image.source="https://github.com/Bio2Byte/Scop3P-notebooks" \
    org.opencontainers.image.licenses="Apache-2.0" \
    org.opencontainers.image.version="${VERSION}" \
    org.opencontainers.image.revision="${VCS_REF}" \
    org.opencontainers.image.created="${BUILD_DATE}"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/apps

WORKDIR /apps

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ 't-coffee' hmmer && \
    rm -rf /var/lib/apt/lists/* && \
    pip install -U pip

COPY requirements-biophysics.txt /app/requirements-biophysics.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install -r /app/requirements-biophysics.txt

COPY requirements-shiny.txt /app/requirements-shiny.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install -r /app/requirements-shiny.txt

COPY --chmod=755 TM-align /usr/local/bin/TM-align
COPY ./apps .

EXPOSE 8000

FROM runtime AS peptide-mapper
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/peptide_mapper/app.py"]

FROM runtime AS structure-viz
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/structure_viz/app.py"]

FROM runtime AS mutation-effect
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/mutation_effect/app.py"]

FROM runtime AS scop3p-toolkit
CMD ["python", "-m", "uvicorn", "portal.main:app", "--host", "0.0.0.0", "--port", "8000"]

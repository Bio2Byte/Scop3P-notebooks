# syntax=docker/dockerfile:1.7
FROM python:3.12-slim AS builder

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION=dev

ENV VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ hmmer t-coffee && \
    python -m venv "$VIRTUAL_ENV" && \
    "$VIRTUAL_ENV/bin/pip" install -U pip && \
    rm -rf /var/lib/apt/lists/*

COPY requirements-biophysics.txt /tmp/requirements-biophysics.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    "$VIRTUAL_ENV/bin/pip" install -r /tmp/requirements-biophysics.txt

COPY requirements-shiny.txt /tmp/requirements-shiny.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    "$VIRTUAL_ENV/bin/pip" install -r /tmp/requirements-shiny.txt

FROM python:3.12-slim AS runtime

LABEL org.opencontainers.image.title="Scop3P-Toolkit" \
    org.opencontainers.image.description="Scop3P-Toolkit: Tools for exploring and extending Scop3P" \
    org.opencontainers.image.url="https://github.com/Bio2Byte/Scop3P-notebooks" \
    org.opencontainers.image.source="https://github.com/Bio2Byte/Scop3P-notebooks" \
    org.opencontainers.image.licenses="Apache-2.0" \
    org.opencontainers.image.version="${VERSION}" \
    org.opencontainers.image.revision="${VCS_REF}" \
    org.opencontainers.image.created="${BUILD_DATE}"

ENV VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/apps

RUN apt-get update && apt-get install -y --no-install-recommends \
    't-coffee' hmmer && \
    rm -rf /var/lib/apt/lists/*

RUN groupadd --system scop3p && \
    useradd --system --gid scop3p --create-home --home-dir /home/scop3p scop3p && \
    mkdir -p /apps /tmp/scop3p && \
    chown -R scop3p:scop3p /apps /tmp/scop3p /home/scop3p

RUN printf '%s\n' 'export VIRTUAL_ENV=/opt/venv' 'export PATH=/opt/venv/bin:$PATH' > /etc/profile.d/scop3p-venv.sh

WORKDIR /apps

COPY --from=builder /opt/venv /opt/venv
COPY --chmod=755 TM-align /usr/local/bin/TM-align
COPY --chown=scop3p:scop3p ./apps .

EXPOSE 8000

USER scop3p

FROM runtime AS peptide-mapper
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/peptide_mapper/app.py"]

FROM runtime AS structure-viz
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/structure_viz/app.py"]

FROM runtime AS mutation-effect
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/mutation_effect/app.py"]

FROM runtime AS scop3p-toolkit
CMD ["python", "-m", "uvicorn", "portal.main:app", "--host", "0.0.0.0", "--port", "8000"]

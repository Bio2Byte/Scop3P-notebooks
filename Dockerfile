# syntax=docker/dockerfile:1
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    VOILA_PORT=8866 \
    VOILA_IP=0.0.0.0

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python deps
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --pre -r /app/requirements.txt

# --- TM-align (optional: copy a known-good binary) ---
# Put your TM-align binary next to this Dockerfile as "TM-align"
COPY TM-align /usr/local/bin/TM-align
RUN chmod +x /usr/local/bin/TM-align || true

# Copy ALL notebooks into the image
COPY notebooks/ /app/notebooks/

EXPOSE 8866

# Serve the directory -> Voila landing page lists notebooks
CMD ["sh", "-c", "\
voila /app/notebooks \
  --Voila.ip=${VOILA_IP} \
  --port=${VOILA_PORT} \
  --no-browser \
  --VoilaConfiguration.enable_nbextensions=True \
  --Voila.tornado_settings='{\"headers\": {\"Content-Security-Policy\": \"frame-ancestors *\"}}' \
"]

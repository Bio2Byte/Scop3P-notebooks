# syntax=docker/dockerfile:1
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
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
RUN pip install --no-cache-dir -r /app/requirements.txt

# --- TM-align (best: copy a known-good binary) ---
# Put your TM-align binary next to this Dockerfile as "TM-align"
# and call it via "TM-align" in your notebook (PATH-based).
# If you don't need TM-align for the teaching flow, you can remove this block.
# Optional include: set build-arg `INCLUDE_TMALIGN=true` when building to
# indicate you are including `TM-align` in the build context. Note: `COPY`
# will fail if the file is missing; remove this block if you don't provide it.
ARG INCLUDE_TMALIGN=false
COPY TM-align /usr/local/bin/TM-align
RUN chmod +x /usr/local/bin/TM-align || true

# App notebook
COPY Scop3P_PTM_structure_viz_voila.ipynb /app/notebook.ipynb

EXPOSE 8866

# Important flags:
# - enable_nbextensions: needed for many widget stacks
# - allow_origin='*' + CSP frame-ancestors: helps PyVis/iframes in some setups
CMD ["sh", "-c", "\
voila /app/notebook.ipynb \
  --Voila.ip=${VOILA_IP} \
  --port=${VOILA_PORT} \
  --base_url=/scop3p-notebook/\
  --no-browser \
  --VoilaConfiguration.enable_nbextensions=True \
  --Voila.tornado_settings='{\"headers\": {\"Content-Security-Policy\": \"frame-ancestors *\"}}' \
"]

# syntax=docker/dockerfile:1.7
#
# Toolkit portal. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make scop3p-toolkit
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: Toolkit portal" \
    org.opencontainers.image.description="All Scop3P-Toolkit apps behind one navbar. The published Galaxy-facing image."

ENV SCOP3P_APP_NAME=scop3p-toolkit

CMD ["python", "-m", "uvicorn", "portal.main:app", "--host", "0.0.0.0", "--port", "8000"]

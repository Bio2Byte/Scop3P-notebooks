# syntax=docker/dockerfile:1.7
#
# RIN Alignment. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make rinalign
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: RIN Alignment" \
    org.opencontainers.image.description="Build, diff and align residue interaction networks from two structures."

ENV SCOP3P_APP_NAME=rinalign

CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/rinalign/app.py"]

# syntax=docker/dockerfile:1.7
#
# Peptide Mapper. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make peptide-mapper
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: Peptide Mapper" \
    org.opencontainers.image.description="Map Scop3P or uploaded phospho-peptides onto AlphaFold structures."

ENV SCOP3P_APP_NAME=peptide-mapper

CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/peptide_mapper/app.py"]

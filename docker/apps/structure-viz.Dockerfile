# syntax=docker/dockerfile:1.7
#
# Structure Visualisation. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make structure-viz
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: Structure Visualisation" \
    org.opencontainers.image.description="PTMs, disease variants, biophysical properties and TM-align on 3D structures."

ENV SCOP3P_APP_NAME=structure-viz

CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/structure_viz/app.py"]

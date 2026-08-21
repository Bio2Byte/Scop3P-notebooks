# syntax=docker/dockerfile:1.7
#
# Topology Viewer. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make topology-viewer
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: Topology Viewer" \
    org.opencontainers.image.description="2D secondary-structure topology diagrams beside a 3D viewer."

ENV SCOP3P_APP_NAME=topology-viewer

CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/topology_viewer/app.py"]

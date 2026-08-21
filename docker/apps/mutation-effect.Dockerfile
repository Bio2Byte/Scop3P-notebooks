# syntax=docker/dockerfile:1.7
#
# Mutation Effect. One app, one image, one CMD -- so which app you get depends on which
# Dockerfile you build, not on stage ordering.
#
# Requires the base image to exist first:
#     make base && make mutation-effect
#
# Build context is this directory; nothing but the base image is needed.
ARG BASE_IMAGE=bio2byte/scop3p-base:local
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Scop3P-Toolkit: Mutation Effect" \
    org.opencontainers.image.description="Compare wild-type and mutant Bio2Byte biophysical predictions."

ENV SCOP3P_APP_NAME=mutation-effect

CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "/apps/mutation_effect/app.py"]

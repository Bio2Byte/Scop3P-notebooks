FROM quay.io/jupyter/minimal-notebook:python-3.12

COPY --chown=root --chmod=755 TM-align /usr/local/bin/TM-align

RUN --mount=type=bind,source=requirements.txt,target=/tmp/requirements.txt pip install -U pip && pip install -r /tmp/requirements.txt

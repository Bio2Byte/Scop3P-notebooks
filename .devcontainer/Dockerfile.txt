FROM mcr.microsoft.com/devcontainers/python:3.11

RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# Install TM-align
RUN wget -q https://zhanggroup.org/TM-align/TMalign && \
    chmod +x TMalign && \
    mv TMalign /usr/local/bin/TMalign && \
    ln -sf /usr/local/bin/TMalign /usr/local/bin/TM-align

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install -U pip && pip install -r /tmp/requirements.txt

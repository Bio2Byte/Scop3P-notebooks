FROM quay.io/jupyter/scipy-notebook:python-3.11

COPY TM-align /usr/local/bin/TM-align
RUN chmod +x /usr/local/bin/TM-align

COPY requirements.txt /tmp/requirements.txt
RUN pip install -U pip && pip install -r /tmp/requirements.txt

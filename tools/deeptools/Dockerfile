FROM python:3.10-slim@sha256:c1e4e6c01eb489c422288b2de34b0761ca316f7a2d98e2c33f47659a73ed108a
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib HOME=/tmp
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.lock /opt/requirements.lock
RUN pip install --no-cache-dir -r /opt/requirements.lock \
    && bamCoverage --version \
    && pip check
RUN apt-get update && apt-get install -y --no-install-recommends samtools \
    && samtools --version > /opt/samtools-version.txt \
    && rm -rf /var/lib/apt/lists/*
COPY worker.py /opt/worker.py
WORKDIR /output
ENTRYPOINT ["python", "/opt/worker.py"]

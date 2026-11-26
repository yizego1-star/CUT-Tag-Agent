FROM cuttag-deeptools:3.5.6
USER root
RUN apt-get update && apt-get install -y --no-install-recommends fastp=0.24.0+dfsg-1 bowtie2=2.5.4-1+b2 bedtools=2.31.1+dfsg-2+b1 macs=3.0.2-1+b3 \
    && dpkg-query -W fastp bowtie2 bedtools macs samtools > /opt/pipeline-versions.txt \
    && fastp --version && bowtie2 --version && /usr/bin/macs3 --version \
    && rm -rf /var/lib/apt/lists/*
COPY worker.py /opt/pipeline-worker.py
ENTRYPOINT ["python", "/opt/pipeline-worker.py"]

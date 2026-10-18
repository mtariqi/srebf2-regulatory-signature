FROM python:3.14-slim

LABEL org.opencontainers.image.title="SREBF2 Phase 12"
LABEL org.opencontainers.image.description="Reproducible SREBF2 regulatory signature statistical analysis"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

WORKDIR /workspace

# Build dependencies for pyBigWig, including remote bigWig support
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        ca-certificates \
        pkg-config && \
    rm -rf /var/lib/apt/lists/*

COPY requirements-phase12.txt /tmp/requirements.txt

RUN python -m pip install --no-cache-dir \
    -r /tmp/requirements.txt && \
    python -m pip check

CMD ["python", "--version"]

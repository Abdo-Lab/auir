# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="MultiQC (v1.0) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl  \
    gcc \
    python-dev \
    && curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py \
    && python /opt/get-pip.py \
    && rm -rf /var/lib/apt/lists/* /opt/get-pip.py

# Setup environment variables
ENV MULTIQC_VERSION="1.0"

# Install MultiQC
RUN pip install multiqc==$MULTIQC_VERSION

# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="BWA (v0.7.15) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV BWA_SOURCE="bwa-0.7.15.tar.bz2" \
    BWA_VERSION="0.7.15" \
    BWA_BIN="bwa" \
    BWA_DEST="/usr/local/bin/bwa"

# Install BWA
RUN wget https://github.com/lh3/bwa/releases/download/v$BWA_VERSION/$BWA_SOURCE -O /opt/$BWA_SOURCE \
    && tar -xvf /opt/$BWA_SOURCE -C /opt \
    && cd /opt/bwa-$BWA_VERSION \
    && make \
    && ln -s /opt/bwa-$BWA_VERSION/$BWA_BIN $BWA_DEST \
    && rm /opt/$BWA_SOURCE

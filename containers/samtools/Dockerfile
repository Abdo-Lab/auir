# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="Samtools (v1.4.1) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV SAMTOOLS_SOURCE="samtools-1.4.1.tar.bz2" \
    SAMTOOLS_VERSION="1.4.1" \
    SAMTOOLS_BIN="samtools" \
    SAMTOOLS_DEST="/usr/local/bin/samtools"

# Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/$SAMTOOLS_SOURCE -O /opt/$SAMTOOLS_SOURCE \
    && tar -xvf /opt/$SAMTOOLS_SOURCE -C /opt \
    && cd /opt/samtools-$SAMTOOLS_VERSION \
    && ./configure \
    && make \
    && ln -s /opt/samtools-$SAMTOOLS_VERSION/$SAMTOOLS_BIN $SAMTOOLS_DEST \
    && rm /opt/$SAMTOOLS_SOURCE

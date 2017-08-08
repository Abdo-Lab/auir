# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="Bedtools2 (v2.26.0) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    python-minimal \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV BEDTOOLS2_SOURCE="bedtools-2.26.0.tar.gz" \
    BEDTOOLS2_VERSION="2.26.0" \
    BEDTOOLS2_BIN="/opt/bedtools2/bin/*" \
    BEDTOOLS2_DEST="/usr/local/bin"

# Install Bedtools2
RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS2_VERSION/$BEDTOOLS2_SOURCE -O /opt/$BEDTOOLS2_SOURCE \
    && tar -xvf /opt/$BEDTOOLS2_SOURCE -C /opt \
    && cd /opt/bedtools2 \
    && make \
    && cp $BEDTOOLS2_BIN $BEDTOOLS2_DEST \
    && rm /opt/$BEDTOOLS2_SOURCE

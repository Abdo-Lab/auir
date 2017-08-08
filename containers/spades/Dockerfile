# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="SPAdes (v3.10.1) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    python \
    wget \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV SPADES_SOURCE="SPAdes-3.10.1-Linux.tar.gz" \
    SPADES_VERSION="3.10.1" \
    SPADES_BIN="bin/spades.py" \
    SPADES_DEST="/usr/local/bin/spades.py"

# Install SPAdes
RUN wget https://github.com/ablab/spades/releases/download/v$SPADES_VERSION/$SPADES_SOURCE -O /opt/$SPADES_SOURCE \
    && tar -xvf /opt/$SPADES_SOURCE -C /opt \
    && cd /opt/SPAdes-$SPADES_VERSION-Linux \
    && ln -s /opt/SPAdes-$SPADES_VERSION-Linux/$SPADES_BIN $SPADES_DEST \
    && rm /opt/$SPADES_SOURCE

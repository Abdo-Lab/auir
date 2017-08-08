# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="QUAST (v4.5) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libfreetype6-dev \
    libpng-dev \
    pkg-config \
    python-matplotlib \
    python-setuptools \
    wget \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV QUAST_SOURCE="quast_4.5.tar.gz" \
    QUAST_VERSION="quast_4.5"

# Install QUAST
RUN wget https://github.com/ablab/quast/archive/$QUAST_SOURCE -O /opt/$QUAST_SOURCE \
    && tar -xvf /opt/$QUAST_SOURCE -C /opt \
    && cd /opt/quast-$QUAST_VERSION \
    && ./setup.py install \
    && rm /opt/$QUAST_SOURCE

# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="Influenza Virus Database Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV IVD_DB="influenza.fna" \
    IVD_DEST="/opt/ivd"

# Download influenza genomes from Influenza Virus Database (IVD)
RUN mkdir $IVD_DEST \
    && wget ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/$IVD_DB -O $IVD_DEST/$IVD_DB

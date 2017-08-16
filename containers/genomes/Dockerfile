FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for influenza virus database genomes" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
 && rm -rf /var/lib/apt/lists/*

ENV IVD_DB="influenza.fna" \
    IVD_DEST="/opt/ivd"

RUN mkdir $IVD_DEST \
    && wget ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/$IVD_DB -O $IVD_DEST/$IVD_DB

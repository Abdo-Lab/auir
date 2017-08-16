FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for QUAST" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

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

ENV QUAST_SOURCE="quast_4.5.tar.gz" \
    QUAST_VERSION="quast_4.5"

RUN wget https://github.com/ablab/quast/archive/$QUAST_SOURCE -O /opt/$QUAST_SOURCE \
    && tar -xvf /opt/$QUAST_SOURCE -C /opt \
    && cd /opt/quast-$QUAST_VERSION \
    && ./setup.py install \
    && rm /opt/$QUAST_SOURCE

FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for Bedtools2" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    python-minimal \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

ENV BEDTOOLS2_SOURCE="bedtools-2.26.0.tar.gz" \
    BEDTOOLS2_VERSION="2.26.0" \
    BEDTOOLS2_BIN="/opt/bedtools2/bin/*" \
    BEDTOOLS2_DEST="/usr/local/bin"

RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS2_VERSION/$BEDTOOLS2_SOURCE -O /opt/$BEDTOOLS2_SOURCE \
    && tar -xvf /opt/$BEDTOOLS2_SOURCE -C /opt \
    && cd /opt/bedtools2 \
    && make \
    && cp $BEDTOOLS2_BIN $BEDTOOLS2_DEST \
    && rm /opt/$BEDTOOLS2_SOURCE

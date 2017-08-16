FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for SPAdes" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    python \
    wget \
 && rm -rf /var/lib/apt/lists/*

ENV SPADES_SOURCE="SPAdes-3.10.1-Linux.tar.gz" \
    SPADES_VERSION="3.10.1" \
    SPADES_BIN="bin/spades.py" \
    SPADES_DEST="/usr/local/bin/spades.py"

RUN wget https://github.com/ablab/spades/releases/download/v$SPADES_VERSION/$SPADES_SOURCE -O /opt/$SPADES_SOURCE \
    && tar -xvf /opt/$SPADES_SOURCE -C /opt \
    && cd /opt/SPAdes-$SPADES_VERSION-Linux \
    && ln -s /opt/SPAdes-$SPADES_VERSION-Linux/$SPADES_BIN $SPADES_DEST \
    && rm /opt/$SPADES_SOURCE

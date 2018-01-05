FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2018-01-05T06:52:21Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for VCFtools" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    dh-autoreconf \
    pkg-config \
    tabix \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

ENV VCF_SOURCE="vcftools-0.1.15.tar.gz" \
    VCF_VERSION="0.1.15" \
    VCF_BIN="src/cpp/vcftools" \
    VCF_DEST="/usr/local/bin/vcftools" \
    PERL5LIB="/opt/vcftools-$VCF_VERSION/src/perl/"

RUN wget https://github.com/vcftools/vcftools/releases/download/v$VCF_VERSION/$VCF_SOURCE -O /opt/$VCF_SOURCE \
    && tar -xvf /opt/$VCF_SOURCE -C /opt \
    && cd /opt/vcftools-$VCF_VERSION \
    && ./autogen.sh \
    && ./configure \
    && make \
    && ln -s /opt/vcftools-$VCF_VERSION/$VCF_BIN $VCF_DEST \
    && rm /opt/$VCF_SOURCE

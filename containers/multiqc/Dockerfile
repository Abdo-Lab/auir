FROM debian:jessie-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for MultiQC" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl  \
    gcc \
    python-dev \
    && curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py \
    && python /opt/get-pip.py \
    && rm -rf /var/lib/apt/lists/* /opt/get-pip.py

ENV MULTIQC_VERSION="1.0"

RUN pip install multiqc==$MULTIQC_VERSION

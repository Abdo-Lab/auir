FROM python:2-slim
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2018-01-07T04:39:25Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for Biopython" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

ENV BIOPYTHON_VERSION="1.70"

RUN pip install biopython==$BIOPYTHON_VERSION

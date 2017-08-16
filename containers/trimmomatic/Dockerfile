FROM openjdk:alpine
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for Trimmomatic" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apk add --update bash && rm -rf /var/cache/apk/*

ENV TRIMMOMATIC_SOURCE="Trimmomatic-0.36.zip" \
    TRIMMOMATIC_HOME="/opt/trimmomatic"

COPY $TRIMMOMATIC_SOURCE $TRIMMOMATIC_HOME/$TRIMMOMATIC_SOURCE

RUN cd $TRIMMOMATIC_HOME \
    && unzip $TRIMMOMATIC_SOURCE \
    && rm $TRIMMOMATIC_SOURCE

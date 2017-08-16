FROM cdeanj/ivd:1.0
MAINTAINER Chris Dean <cdean11@colostate.edu>

LABEL org.label-schema.build-date="2017-08-16T01:31:38Z" \
      org.label-schema.name="AUIR" \
      org.label-schema.description="A Docker image for Blast" \
      org.label-schema.url="https://github.com/cdeanj/auir" \
      org.label-schema.vcs-ref="71f505cbdd9a4501ecbc2a39b2352ddca7b37fb7" \
      org.label-schema.vcs-url="https://github.com/cdeanj/auir" \
      org.label-schema.schema-version="1.0.0-rc1" \
      org.label-schema.vendor="Abdo Lab"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

ENV BLASTDB_DIR="BLASTDB" \
    BLASTDB_TYPE="nucl" \
    BLASTDB_NAME="InfluenzaDB" \
    BLASTDB_INPUT="influenza.fna"

ENV IVD_DIR="ivd"

RUN cd /opt \
    && mkdir $BLASTDB_DIR

ENV BLASTDB "$BLASTDB:/opt/BLASTDB"

RUN cd /opt/$IVD_DIR \
    && makeblastdb -in $BLASTDB_INPUT -dbtype $BLASTDB_TYPE -out $BLASTDB_NAME \
    && mv $BLASTDB_NAME* /opt/$BLASTDB_DIR

# Set the base image
FROM cdeanj/ivd:1.0

# Set image metadata
LABEL author="Christopher Dean" \
      description="NCBI-BLAST+ Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV BLASTDB_DIR="BLASTDB" \
    BLASTDB_TYPE="nucl" \
    BLASTDB_NAME="InfluenzaDB" \
    BLASTDB_INPUT="influenza.fna"

ENV IVD_DIR="ivd"

# Make blast database
RUN cd /opt \
    && mkdir $BLASTDB_DIR

ENV BLASTDB "$BLASTDB:/opt/BLASTDB"

RUN cd /opt/$IVD_DIR \
    && makeblastdb -in $BLASTDB_INPUT -dbtype $BLASTDB_TYPE -out $BLASTDB_NAME \
    && mv $BLASTDB_NAME* /opt/$BLASTDB_DIR

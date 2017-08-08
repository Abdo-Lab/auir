# Set the base image
FROM openjdk:alpine

# Set image metadata
LABEL author="Christopher Dean" \
      description="Trimmomatic (v0.36) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core packages
RUN apk add --update bash && rm -rf /var/cache/apk/*

# Setup environment variables
ENV TRIMMOMATIC_SOURCE="Trimmomatic-0.36.zip" \
    TRIMMOMATIC_HOME="/opt/trimmomatic"

# Add Trimmomatic source file to build context
COPY $TRIMMOMATIC_SOURCE $TRIMMOMATIC_HOME/$TRIMMOMATIC_SOURCE

# Install Trimmomatic
RUN cd $TRIMMOMATIC_HOME \
    && unzip $TRIMMOMATIC_SOURCE \
    && rm $TRIMMOMATIC_SOURCE

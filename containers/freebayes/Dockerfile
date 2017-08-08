# Set the base image
FROM linuxbrew/linuxbrew

# Set image metadata
LABEL author="Christopher Dean" \
      description="FreeBayes (v1.1.0) Docker image" \
      maintainer="cdean11@colostate.edu"

RUN brew tap homebrew/science \
    && brew install freebayes

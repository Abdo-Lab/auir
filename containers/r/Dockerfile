# Set the base image
FROM r-base

# Set image metadata
LABEL author="Christopher Dean" \
      description="R Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN echo 'install.packages(c("data.table", "reshape2", "ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
    && Rscript /tmp/packages.R

# Base image with R
FROM rocker/r-ver:4.3.1

# Install system dependencies for DESeq2 and ggplot2
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libpng-dev \
    libxt-dev \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor manager and packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('DESeq2', 'ggplot2'), ask=FALSE, update=TRUE)"

# Base image with R
FROM rocker/r-ver:4.3.1

# Install system dependencies for R packages + Python
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libpng-dev \
    libxt-dev \
    libz-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor manager and R packages including IRkernel
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('DESeq2', 'ggplot2', 'IRkernel'), ask=FALSE, update=TRUE)"

# Install Python Jupyter + extensions
RUN pip3 install --no-cache-dir \
    jupyter \
    notebook \
    jupyterlab \
    ipykernel \
    jupyter_contrib_nbextensions \
    jupyterlab-git \
    jupyterlab_widgets \
    matplotlib \
    pandas \
    numpy \
    scipy \
    scikit-learn

# Register R kernel with Jupyter (system-wide)
RUN R -e "IRkernel::installspec(user = FALSE)"

# Use an R base image
FROM rocker/r-ver:4.3.1

# Install required system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    libglpk-dev \
    bash \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor package manager
RUN R -e "install.packages('BiocManager')"

# Install DESeq2, clusterProfiler, and Seurat (from CRAN)
RUN R -e "install.packages('Seurat', version = '5.0.1')"
RUN R -e "install.packages('SeuratObject', version = '5.0.1')"
RUN R -e "install.packages('Matrix', version = '1.6-3')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('pheatmap')"

# Optionally, clean up to reduce image size
RUN rm -rf /tmp/*

# Set the working directory
WORKDIR /root

# Start from rocker RStudio image
FROM rocker/rstudio:4.4.1

# Install system dependencies (for compiling R packages)
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pak package manager
RUN R -e "install.packages(c('metrics','svglite', 'BiocManager','funkyheatmap','tidyverse', 'dplyr', 'tidyr', 'stringr','ggplot2', 'ggpubr', 'patchwork', 'pROC', 'devtools', 'knitr', 'rmarkdown', 'remotes'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'circlize'))"

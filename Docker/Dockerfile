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
RUN R -e "install.packages(c('pak','metrics'))"

# Install CRAN packages
RUN R -e "pak::pkg_install(c('tidyverse', 'dplyr', 'tidyr', 'stringr','ggplot2', 'ggpubr', 'patchwork', 'pROC', 'devtools', 'knitr', 'rmarkdown', 'remotes', 'funkyheatmap', 'svglite', 'BiocManager'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'circlize'))"

# Start from rocker RStudio image
FROM rocker/rstudio:4.4.1

# Install system dependencies (for compiling R and Bioconductor packages)
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c('svglite', 'BiocManager', 'funkyheatmap', 'tidyverse', 'dplyr', 'tidyr', 'stringr', 'ggplot2', 'ggpubr', 'patchwork', 'pROC', 'devtools', 'knitr', 'rmarkdown', 'remotes', 'tidytext', 'car', 'PMCMRplus', 'tibble', 'purrr', 'ggrepel', 'RColorBrewer', 'circlize'))"

# Install Bioconductor packages (release matching R 4.4)
RUN R -e "BiocManager::install(version = '3.20', ask = FALSE, update = FALSE); BiocManager::install(c('ComplexHeatmap', 'RnBeads'), ask = FALSE, update = FALSE)"

# Install MeDeCom from GitHub (used by resources/analyze_*_lcms.R)
RUN R -e "remotes::install_github('lutsik/MeDeCom', upgrade = 'never')"

# Fail the build if any required package is missing
# (install.packages only warns when an installation fails)
RUN R -e "pkgs <- c('svglite', 'funkyheatmap', 'tidyverse', 'dplyr', 'tidyr', 'stringr', 'ggplot2', 'ggpubr', 'patchwork', 'pROC', 'knitr', 'rmarkdown', 'tidytext', 'car', 'PMCMRplus', 'tibble', 'purrr', 'ggrepel', 'RColorBrewer', 'circlize', 'ComplexHeatmap', 'RnBeads', 'MeDeCom'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop('Missing packages: ', paste(missing, collapse = ', '))"

# Benchmarking Study – Reproducibility Repository

This repository contains the **R Markdown notebooks** (and corresponding **HTML outputs**) used to reproduce the analyses from our benchmarking study [publication link to be added].

The tutorials included here guide you through the main analyses presented in the paper. For additional visualizations published in the **Supplementary Figures** or made available in the **R Shiny App**, please visit:

* 📊 [Interactive App](https://sunny.cmb.ugent.be/3fy5CTR4gXjcKMHj0zz1bxsGEEkHVsOnC8BjXYT5miRB0QGwid/)
* 💻 [App Source Code](https://github.com/VIBTOBIlab/Shiny-DNAmBenchmarking.git)

---

## Repository Structure

```
.
├── Docker/       # Contains Dockerfile used to build the Docker Container
├── notebooks/    # R Markdown (.Rmd) files
│   └── html/     # Rendered HTML files from the notebooks
├── R/            # Custom R functions used across the notebooks
├── resources/    # Results, datasets, and supporting files for analyses & plots
└── plots/        # Generated automatically by the notebooks, organized by step
```

---

## Run the analyses using a Docker container

You can run the full analysis inside a Docker container (already built using the following [Dockerfile](Docker/Dockerfile)) without manually installing dependencies:

1. Run the container, mounting the current repository:

   ```bash
   docker run --rm \                                    
      -p 8787:8787 \
      -e PASSWORD=mypassword \ # Modify the password with a personal one
      -v $(pwd):/home/rstudio/benchmark \
      egiuili/benchmark-rstudio:v1
   ```
2. Open http://localhost:8787 in your web browser to access RStudio Server inside the container.

   - Default username: rstudio
   - Default password: mypassword

Once logged in, you can open and run the .Rmd notebooks located in /home/rstudio/benchmark/notebooks/.

---

## Run the analyses with locally installed dependencies

Alternatively, you can run the analysis using your local environment, but make sure that the following R packages are installed before running the notebooks.

```r
# Install CRAN packages
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "stringr",
  "ggplot2", "ggpubr", "patchwork", "pROC",
  "metrics", "devtools", "knitr", "rmarkdown", 
  "remotes", "funkyheatmap", "svglite"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "circlize"))
```

---

## Analysis Workflow

The analyses should be run in the following order. Each notebook generates plots saved in a dedicated subfolder under `plots/`:

1. **01_heatmaps_all_combinations.Rmd**

   * Generates heatmaps of all possible in silico mixture combinations.
   * Outputs to `plots/01_comprehensive_heatmaps/`.

2. **02_limit_of_detection.Rmd**

   * Reproduces the **limit of detection (LoD)** results.
   * Outputs to `plots/02_lod_heatmaps/`.

3. **03_sequencing_depth_robustness.Rmd**

   * Evaluates **AUC-ROC robustness** across sequencing depths.
   * Outputs to `plots/03_seqdepth_robustness/`.

4. **04_reference_free_analysis.Rmd**

   * Compares **reference-based vs. reference-free** tools.
   * Outputs to `plots/04_refree_vs_refbased/`.

5. **05_additional_plots.Rmd**

   * Generates additional visualizations:

     * AUC-ROC tool rankings (DotPlots)
     * Ranking distributions
     * Dataset-specific rankings
   * Outputs to `plots/05_additional_plots/`.

6. **06_funkyheatmap.Rmd**

   * Produces the **funkyheatmap** figure from the publication.
   * Outputs to `plots/06_funkyheatmap/`.

7. **07_scalability.Rmd**

   * Analyzes and visualizes **scalability results**.
   * Outputs to `plots/07_scalability/`.

---

## Summary of Output Structure

After running the notebooks, the `plots/` folder will be organized as follows:

```
plots/
├── 01_comprehensive_heatmaps/   # Heatmaps of all in silico mixtures
├── 02_lod_heatmaps/             # Limit of detection plots
├── 03_seqdepth_robustness/      # Sequencing depth robustness plots
├── 04_refree_vs_refbased/       # Reference-free vs. reference-based comparison plots
├── 05_additional_plots/         # Ranking and additional visualizations
├── 06_funkyheatmap/             # Funkyheatmap figure from publication
└── 07_scalability/              # Scalability analysis plots
```

---

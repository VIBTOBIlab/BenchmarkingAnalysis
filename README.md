# Benchmarking Study – Reproducibility Repository

This repository contains the **R Markdown notebooks** (and corresponding **HTML outputs**) used to reproduce the analyses from our benchmarking study [publication link to be added].

The tutorials included here guide you through the main analyses presented in the paper. For additional visualizations published in the **Supplementary Figures** or made available in the **R Shiny App**, please visit:

* 📊 [Interactive App](https://sunny.cmb.ugent.be/3fy5CTR4gXjcKMHj0zz1bxsGEEkHVsOnC8BjXYT5miRB0QGwid/)
* 💻 [App Source Code](https://github.com/VIBTOBIlab/Shiny-DNAmBenchmarking.git)

---

## 📚 Table of Contents

1. [Repository Structure](#repository-structure)
2. [Run the R markdown notebooks](#run-the-r-markdown-notebooks)
   - [Run the analyses using Docker](#run-the-analyses-using-a-docker-container)
   - [Run the analyses locally](#run-the-analyses-with-locally-installed-dependencies)
   - [Analysis workflow](#analysis-workflow)
   - [Summary of output structure](#summary-of-output-structure)
3. [Calculate the sequencing saturation curves](#sequencing-saturation-calculation)
   - [What is a sequencing saturation curve?](#sequencing-saturation-calculation)
   - [Run the sequencing saturation calculation](#31-how-to-run-a-sequencing-saturation-calculation)


---

## 1. Repository Structure

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
## 2. Run the R markdown notebooks

### 2.1 Run the analyses using a Docker container

You can run the full analysis inside a Docker container (already built using the following [Dockerfile](Docker/Dockerfile)) without manually installing dependencies:

1. Run the container, mounting the current repository:

   ```bash
   docker run --rm \                                    
      -p 8787:8787 \
      -e PASSWORD=mypassword \ # Modify the password with a personal one
      -v "$(pwd)":/home/rstudio/benchmark \
      egiuili/benchmark-rstudio:v3
   ```

2. Open [http://localhost:8787](http://localhost:8787) in your web browser to access RStudio Server inside the container.

   * Default username: `rstudio`
   * Default password: `mypassword`

Once logged in, you can open and run the `.Rmd` notebooks located in `/home/rstudio/benchmark/notebooks/`.

---

### 2.2 Run the analyses with locally installed dependencies

Alternatively, you can run the analysis using your local environment, but make sure that the following R packages are installed before running the notebooks.

```r
# Install CRAN packages
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "stringr",
  "ggplot2", "ggpubr", "patchwork", "pROC",
  "metrics", "devtools", "knitr", "rmarkdown", 
  "remotes", "funkyheatmap", "svglite", "tidytext"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "circlize"))
```

---

### 2.3 Analysis Workflow

The analyses should be run in the following order. Each notebook generates plots saved in a dedicated subfolder under `plots/`:

1. **01_heatmaps_all_combinations.Rmd**

   * Generates heatmaps of all possible *in silico* mixture combinations.
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

8. **08_preciseness_tumor_fractions.Rmd**

   * Analyzes and visualizes **RMSE and SCC for tumor fraction prediction accuracy**.
   * Outputs to `plots/08_preciseness/`.

---

### 2.4 Summary of Output Structure

After running the notebooks, the `plots/` folder will be organized as follows:

```
plots/
├── 01_comprehensive_heatmaps/   # Heatmaps of all in silico mixtures
├── 02_lod_heatmaps/             # Limit of detection plots
├── 03_seqdepth_robustness/      # Sequencing depth robustness plots
├── 04_refree_vs_refbased/       # Reference-free vs. reference-based comparison plots
├── 05_additional_plots/         # Ranking and additional visualizations
├── 06_funkyheatmap/             # Funkyheatmap figure from publication
├── 07_scalability/              # Scalability analysis plots
└── 08_preciseness/              # Preciseness analysis plots
```

---

## 3. What is a sequencing saturation curve?

To calculate the **sequencing saturation** of an RRBS sample, we adopted the following strategy:

For each unique sample (e.g., an *in silico* mixture composed of 10% tumor reads and 90% healthy reads), we computed the number of **unique CpGs covered by at least 3 reads** at five different sequencing depths: **2M, 5M, 10M, 15M, and 20M**.

We then fit the following curve using the `scipy.optimize.curve_fit` function:

$$
y = \beta_0 \cdot \arctan(\beta_1 \cdot x)
$$

We chose the **arctangent function** because it exhibits an **asymptotic growth** similar to sequencing saturation.
For large values of $x$ (as $x \to \infty$), the asymptote can be computed as:

$$
\text{asymptote} = \beta_0 \cdot \frac{\pi}{2}
$$

The **sequencing saturation value** at each depth was then calculated as:

$$
\text{Saturation} = \frac{\text{Number of unique CpGs (≥3 counts)}}{\text{Asymptote}}
$$

This approach allows estimation of the theoretical **maximum number of CpGs** that can be detected given an infinite sequencing depth, and quantifies how close each sample is to reaching sequencing saturation.

---

### 3.1 How to run a sequencing saturation calculation

You can reproduce the sequencing saturation analysis using the provided Python script:

```bash
# If not, install matplotlib
pip install numpy pandas matplotlib scipy
# Run the script
python3 resources/plot_seqstat.py \
    --cpgs_file resources/cov_cpgcounts_rrbs.csv \
    --read_file resources/bam_readcounts_rrbs.csv \
    --percentages 0.1,0.25,0.50,0.75,1
```
Where the [**--cpgs_file**](resources/cov_cpgcounts.csv) corresponds to a csv file with sample name, percentage of downsampling, minimum number of counts per CpGs used and number of CpGs found the the corresponding COV file. To calculate the number of CpGs in a COV file you can run the following script:

```bash
zcat <your_cov.gz_file> | awk -v OFS='\t' -v i="$i" '$5 + $6 >= 3' | wc -l
```
Where 3 corresponds to the minimum number of counts per CpGs.

The [**--read_file**](resources/bam_readcounts.csv) corresponds to a csv file with sample name, percentage of downsampling and number of reads present in that BAM samples.
#!/usr/bin/env Rscript
# ===============================================================
# Compare MeDeCom latent methylation components (LMCs) to reference profiles
#
# Loads a MeDeCom result (MeDeComSet, as returned by runMeDeCom) and a
# reference matrix of known cell-type methylation profiles (e.g. healthy and
# tumor), then plots how well the inferred LMCs match the references
# (dendrogram + heatmap) and the estimated proportions per sample.
#
# Usage (from the repository root):
#   Rscript resources/analyze_medecom_lcms.R \
#       [--medecom=resources/medecom_lcm.rds] \
#       [--reference=resources/reference_samples_for_lcm.csv] \
#       [--K=2] [--lambda=0] \
#       [--out=plots/lcm_analysis/medecom_lcms.pdf]
#
# Author: Edoardo Giuili
# ===============================================================

suppressPackageStartupMessages(library(MeDeCom))

# ---------------------------------------------------------------
# Command line arguments (--key=value), with defaults
# ---------------------------------------------------------------
parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  for (a in args) {
    kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(kv) != 3 || !kv[2] %in% names(defaults)) {
      stop("Unknown argument: ", a, "\nValid arguments: ",
           paste0("--", names(defaults), collapse = ", "))
    }
    defaults[[kv[2]]] <- kv[3]
  }
  defaults
}

opt <- parse_args(list(
  medecom   = "resources/medecom_lcm.rds",
  reference = "resources/reference_samples_for_lcm.csv",
  K         = "2",
  lambda    = "0",
  out       = "plots/lcm_analysis/medecom_lcms.pdf"
))
K <- as.numeric(opt$K)
lambda <- as.numeric(opt$lambda)

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------

# Impute missing values with the mean of their row (region);
# rows that are entirely NA fall back to the global mean
impute_row_means <- function(m) {
  row_means <- rowMeans(m, na.rm = TRUE)
  idx <- which(is.na(m), arr.ind = TRUE)
  m[idx] <- row_means[idx[, 1]]
  m[is.na(m)] <- mean(m, na.rm = TRUE)
  m
}

# ---------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------
obj <- readRDS(opt$medecom)
cat("MeDeCom result:", obj@dataset_info$m, "regions x",
    obj@dataset_info$n, "samples\n")
cat("Available K:", paste(obj@parameters$Ks, collapse = ", "),
    "| lambdas:", paste(obj@parameters$lambdas, collapse = ", "), "\n")
if (!K %in% obj@parameters$Ks) stop("K = ", K, " not found in the MeDeCom result")
if (!lambda %in% obj@parameters$lambdas) stop("lambda = ", lambda, " not found in the MeDeCom result")

# Reference: regions x cell types, first column holds the region IDs
ref <- read.csv(opt$reference, row.names = 1, check.names = FALSE)
ref <- as.matrix(ref)
storage.mode(ref) <- "numeric"
ref <- impute_row_means(ref)

# ---------------------------------------------------------------
# Align reference and LMCs on the same regions
# ---------------------------------------------------------------
# plotLMCs() compares LMCs and references row by row, so both matrices must
# describe the same regions in the same order
lmcs <- getLMCs(obj, K = K, lambda = lambda)
if (!is.null(rownames(lmcs))) {
  common <- intersect(rownames(lmcs), rownames(ref))
  cat(length(common), "of", nrow(lmcs), "regions found in the reference\n")
  if (length(common) == 0) stop("No shared regions between MeDeCom LMCs and the reference")
  # Restrict the LMCs of every (K, lambda) of the first CpG subset to the shared regions
  T_cells <- obj@outputs[[1]]$T
  for (i in seq_along(T_cells)) T_cells[[i]] <- T_cells[[i]][common, , drop = FALSE]
  obj@outputs[[1]]$T <- T_cells
  obj@dataset_info$m <- length(common)
  ref <- ref[common, , drop = FALSE]
} else if (nrow(lmcs) != nrow(ref)) {
  # Without region names the rows cannot be matched safely
  stop("The MeDeCom LMCs have no region names and ", nrow(lmcs),
       " rows, while the reference has ", nrow(ref), " rows.\n",
       "Provide a reference with the same regions (in the same order) ",
       "as the MeDeCom input, or add region names to the LMCs.")
}

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
pdf(opt$out, width = 8, height = 7)

# LMCs vs reference profiles
plotLMCs(obj, K = K, lambda = lambda, type = "dendrogram", Tref = ref, center = TRUE)
plotLMCs(obj, K = K, lambda = lambda, type = "heatmap", Tref = ref)

# Estimated proportions of each LMC per sample
plotProportions(obj, type = "barplot", K = K, lambda = lambda)
plotProportions(obj, type = "heatmap", K = K, lambda = lambda)

invisible(dev.off())
cat("Plots saved to", opt$out, "\n")

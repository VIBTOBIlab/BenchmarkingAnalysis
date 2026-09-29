#!/usr/bin/env Rscript
# ===============================================================
# Compare RefFreeCellMix latent methylation components (LMCs) to reference profiles
#
# RefFreeCellMix returns two matrices:
#   - W (Mu): regions x K methylation profiles of the latent components
#   - H (Omega): samples x K estimated proportions
# This script wraps them into a MeDeComSet, so the MeDeCom plotting functions
# can be reused, then plots how well the LMCs match known reference profiles
# (dendrogram + heatmap) and the estimated proportions per sample.
#
# Usage (from the repository root):
#   Rscript resources/analyze_refreecellmix_lcms.R \
#       [--lcms=resources/refreecellmix_lcm.csv] \
#       [--proportions=resources/refreecellmix_estimated_proportions.csv] \
#       [--reference=resources/reference_samples_for_lcm.csv] \
#       [--out=plots/lcm_analysis/refreecellmix_lcms.pdf]
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
  lcms        = "resources/refreecellmix_lcm.csv",
  proportions = "resources/refreecellmix_estimated_proportions.csv",
  reference   = "resources/reference_samples_for_lcm.csv",
  out         = "plots/lcm_analysis/refreecellmix_lcms.pdf"
))

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

read_matrix <- function(path) {
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "numeric"
  m
}

# ---------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------
W   <- read_matrix(opt$lcms)         # regions x K
H   <- read_matrix(opt$proportions)  # samples x K
ref <- impute_row_means(read_matrix(opt$reference))  # regions x cell types

K <- ncol(W)
lambda <- 0  # RefFreeCellMix has no lambda, 0 is used as placeholder
if (ncol(H) != K) stop("LMCs (", K, ") and proportions (", ncol(H), ") have a different number of components")
colnames(W) <- colnames(H) <- paste0("LMC", seq_len(K))

# ---------------------------------------------------------------
# Align reference and LMCs on the same regions
# ---------------------------------------------------------------
# plotLMCs() compares LMCs and references row by row, so both matrices must
# describe the same regions in the same order
common <- intersect(rownames(W), rownames(ref))
cat(length(common), "of", nrow(W), "LMC regions found in the reference (",
    nrow(ref), "reference regions )\n")
if (length(common) == 0) stop("No shared regions between the LMCs and the reference")
W <- W[common, , drop = FALSE]
ref <- ref[common, , drop = FALSE]

# ---------------------------------------------------------------
# Build a MeDeComSet with the same structure as runMeDeCom() output
# ---------------------------------------------------------------
# outputs are indexed by CpG subset (a single one here, "1"); T and A are
# list-matrices indexed by (K, lambda) holding:
#   T: regions x K (W), A: K x samples (t(H))
as_cell <- function(x) {
  matrix(list(x), 1, 1, dimnames = list(paste0("K_", K), paste0("lambda_", lambda)))
}
obj <- MeDeComSet(
  parameters   = list(cg_subsets = c("All CpGs" = 1L), Ks = K, lambdas = lambda),
  outputs      = list(`1` = list(T = as_cell(W), A = as_cell(t(H)))),
  dataset_info = list(m = nrow(W), n = nrow(H))
)

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

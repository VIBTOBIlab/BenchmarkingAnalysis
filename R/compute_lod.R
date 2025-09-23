compute_lod <- function(data, comparisons, correction.method = "BH") {
  #-------------------------------------------------------------------
  # Perform pairwise Wilcoxon rank-sum tests with multiple test correction
  #
  # Args:
  #   data: dataframe containing at least two columns:
  #         - expected_tf: grouping variable (factor/character)
  #         - predicted_tf: numeric values to compare
  #   comparisons: list of length-2 character vectors specifying groups to compare
  #   correction.method: string, method used for p.adjust (e.g. "BH", "bonferroni")
  #
  # Returns:
  #   A dataframe with:
  #     - group1, group2: compared groups
  #     - p: raw p-value
  #     - p.adj: adjusted p-value
  #     - y.position: suggested y-position for annotation
  #     - star: significance stars
  #-------------------------------------------------------------------
  
  results <- list()
  
  # --- Run Wilcoxon tests for each specified comparison ---
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Subset only rows for the two groups
    subset_data <- data[data$expected_tf %in% comp, ]
    
    # Perform one-sided Wilcoxon rank-sum test
    test <- wilcox.test(predicted_tf ~ expected_tf,
                        data = subset_data,
                        alternative = "less")
    
    # Collect results
    results[[length(results) + 1]] <- data.frame(
      group1 = group1,
      group2 = group2,
      p = test$p.value
    )
  }
  
  # Combine results into a single dataframe
  stats <- do.call(rbind, results)
  
  # --- Multiple testing correction ---
  stats$p.adj <- p.adjust(stats$p, method = correction.method)
  
  # --- Compute y-axis annotation positions ---
  num_comparisons <- nrow(stats)
  max_val <- max(data$predicted_tf, na.rm = TRUE)
  
  stats <- stats %>%
    arrange(as.numeric(group2)) %>%
    mutate(
      y.position = seq(from = max_val * 1.1,      # start slightly above max
                       to   = max_val * 1.8,      # spread comparisons vertically
                       length.out = num_comparisons),
      star = case_when(
        p.adj <= 0.0001 ~ "****",
        p.adj <= 0.001  ~ "***",
        p.adj <= 0.01   ~ "**",
        TRUE            ~ "ns"   # or "" to hide non-significant
      )
    )
  
  # --- Format p-values for readability ---
  stats <- stats %>%
    mutate(
      p      = formatC(p, format = "e", digits = 1),
      p.adj  = formatC(p.adj, format = "e", digits = 1)
    )
  
  return(stats)
}

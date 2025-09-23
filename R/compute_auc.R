# Function to compute the ROC curve and AUC for predicted vs expected tumor fractions
# Arguments:
#   filt_df  : A data frame containing at least two columns:
#                - expected_tf: the expected tumor fraction (ground truth)
#                - predicted_tf: the predicted tumor fraction (output of deconvolution)
#   fraction : The tumor fraction corresponding to the current evaluation (numeric or factor)
#   deconv   : Name of the deconvolution tool used (character)
# Returns:
#   A data frame with the ROC curve (FPR, TPR, thresholds) and AUC, along with the 
#   fraction and tool name for annotation purposes

roc.obj <- function(true_labels, predicted_scores) {
  
  # Transform all the expected fractions > 0 to class 1
  true_labels[true_labels>0] <- 1
  
  # Compute the AUC-ROC 
  roc_obj <- roc(true_labels, predicted_scores, quiet = TRUE)
  
  return(roc_obj)
}

compute_auc <- function(filt_df, fraction, deconv) {
  
  # Compute the ROC curve object using expected vs predicted tumor fractions
  roc_curve <- roc.obj(filt_df$expected_tf, filt_df$predicted_tf)
  
  # Construct a data frame with ROC metrics and additional metadata
  tmp <- data.frame(
    fpr = 1 - rev(roc_curve$specificities),      # False Positive Rate (1 - specificity), reversed to match increasing thresholds
    tpr = rev(roc_curve$sensitivities),          # True Positive Rate (sensitivity), reversed for increasing thresholds
    thresholds = rev(roc_curve$thresholds),      # Thresholds used to compute TPR and FPR, reversed
    auc = rev(roc_curve$auc),                    # Area Under the Curve, repeated/reversed to match thresholds
    expected_tf = fraction,                      # Tumor fraction for this evaluation
    deconv_tool = deconv,                        # Name of the deconvolution tool
    dmr_tool = unique(filt_df$dmr_tool),         # Name of the DMR tool
    seq_depth = unique(filt_df$seq_depth),       # Sequencing depth analysed
    tumor_type = unique(filt_df$tumor_type)      # Tumor type analysed
  )
  
  # Return the annotated ROC curve and AUC data frame
  return(tmp)
}


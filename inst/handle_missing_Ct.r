# clean <- handle_missing_Ct(data_2factor3ref,
#                            target_ct.cols = c( "Ct_PO","Ct_GAPDH", "Ct_ref2"),
#                            ref_ct.cols = c("Ct_ref3"),
#                            factor.cols = c("factor1", "factor2"),
#                            replicate.col = "Rep",
#                            method = "mean",
#                            fixed.value = 40)
# clean


handle_missing_Ct <- function(
    x,
    target_ct.cols = NULL,
    ref_ct.cols = NULL,
    factor.cols,
    replicate.col,
    method = c("mean", "fixed", "drop"),
    fixed.value = 40
) {
  
  method <- match.arg(method)
  
  #checks
  if (!is.data.frame(x))
    stop("`x` must be a data.frame")
  
  all_ct_cols <- c(target_ct.cols, ref_ct.cols)
  if (!all(all_ct_cols %in% colnames(x)))
    stop("Some Ct columns are missing")
  
  if (!all(factor.cols %in% colnames(x)))
    stop("Some factor columns are missing")
  
  if (!replicate.col %in% colnames(x))
    stop("Replicate column not found")
  
  # DROP
  if (method == "drop") {
    keep <- complete.cases(x[, all_ct_cols, drop = FALSE])
    return(x[keep, , drop = FALSE])
  }
  
  # FIXED
  if (method == "fixed") {
    # Target genes get fixed value
    if (!is.null(target_ct.cols)) {
      for (ct in target_ct.cols) {
        x[[ct]][is.na(x[[ct]])] <- fixed.value
      }
    }
    # Reference genes also optionally fixed (or keep NA)
    if (!is.null(ref_ct.cols)) {
      for (ct in ref_ct.cols) {
        x[[ct]][is.na(x[[ct]])] <- fixed.value
      }
    }
    return(x)
  }
  
  # MEAN (for reference genes only, by factor-level combinations)
  if (method == "mean") {
    
    if (!is.null(target_ct.cols)) {
      # Target genes get fixed value
      for (ct in target_ct.cols) {
        x[[ct]][is.na(x[[ct]])] <- fixed.value
      }
    }
    
    if (!is.null(ref_ct.cols)) {
      # factor-level combinations ONLY (replicate NOT included)
      group_id <- interaction(x[, factor.cols, drop = FALSE], drop = TRUE)
      
      for (ct in ref_ct.cols) {
        
        na_idx <- which(is.na(x[[ct]]))
        if (length(na_idx) == 0) next
        
        for (i in na_idx) {
          
          g <- group_id[i]
          
          # explicitly exclude the current replicate
          vals <- x[[ct]][group_id == g & x[[replicate.col]] != x[[replicate.col]][i]]
          vals <- vals[!is.na(vals)]
          
          if (length(vals) > 0) {
            x[[ct]][i] <- mean(vals)
          } else {
            x[[ct]][i] <- NA
          }
        }
      }
    }
    
    return(x)
  }
}

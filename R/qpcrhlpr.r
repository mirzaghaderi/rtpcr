.convert_to_character <- function(numbers) {
  characters <- rep(" ", length(numbers))
  ok <- !is.na(numbers)
  characters[ok & numbers >= 0.001 & numbers < 0.01] <- "**"
  characters[ok & numbers >= 0.01  & numbers < 0.05] <- "*"
  characters[ok & numbers >= 0.05  & numbers < 0.1 ] <- "."
  return(characters)
}




.cleanup <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    block = NULL,
    set_missing_target_Ct_to_40 = FALSE  # Default FALSE = NA, TRUE = 40
) {
  
  x <- as.data.frame(x)
  
  # Set missing value based on parameter
  missing_value <- if (set_missing_target_Ct_to_40) 40 else NA
  
  # Number of non-gene columns
  n_meta <- numOfFactors + if (!is.null(block)) 2 else 1
  
  # Total gene columns
  gene_cols <- (n_meta + 1):ncol(x)
  
  if (length(gene_cols) %% 2 != 0) stop("Gene columns must be in E/Ct pairs")
  
  # Split into pairs
  gene_pairs <- split(gene_cols, ceiling(seq_along(gene_cols) / 2))
  
  # Reference vs target genes
  ref_pairs <- tail(gene_pairs, numberOfrefGenes)
  target_pairs <- utils::head(gene_pairs, length(gene_pairs) - numberOfrefGenes)
  
  # Helper: check invalid (NA or any string/non-numeric)
  is_invalid <- function(v) {
    is.na(v) | is.na(suppressWarnings(as.numeric(v)))
  }
  
  # Reference genes: NA for invalid values 
  for (pair in ref_pairs) {
    E_col <- pair[1]
    Ct_col <- pair[2]
    
    bad_E <- is_invalid(x[[E_col]])
    x[[E_col]][bad_E] <- NA
    x[[E_col]] <- as.numeric(x[[E_col]])
    
    bad_Ct <- is_invalid(x[[Ct_col]])
    x[[Ct_col]][bad_Ct] <- NA  # Always NA for reference genes
    x[[Ct_col]] <- as.numeric(x[[Ct_col]])
  }
  
  # Target genes - use parameter choice ----
  for (pair in target_pairs) {
    E_col  <- pair[1]
    Ct_col <- pair[2]
    
    bad_E <- is_invalid(x[[E_col]])
    x[[E_col]][bad_E] <- NA
    x[[E_col]] <- as.numeric(x[[E_col]])
    
    bad_Ct <- is_invalid(x[[Ct_col]])
    x[[Ct_col]][bad_Ct] <- missing_value  # NA or 40 based on parameter
    x[[Ct_col]] <- as.numeric(x[[Ct_col]])
  }
  
  return(x)
}






.wide_to_long <- function(df) {
  
  if (ncol(df) < 6) {
    stop("Data frame must contain at least 6 columns.")
  }
  
  # metadata (first two columns)
  meta <- df[, 1:2, drop = FALSE]
  
  # remaining columns (paired by position)
  data_cols <- df[, -(1:2), drop = FALSE]
  
  if (ncol(data_cols) %% 2 != 0) {
    stop("After the first two columns, remaining columns must be in pairs.")
  }
  
  n_pairs <- ncol(data_cols) / 2
  
  out <- do.call(
    rbind,
    lapply(seq_len(n_pairs), function(i) {
      
      e_col  <- data_cols[, 2*i - 1]
      ct_col <- data_cols[, 2*i]
      
      gene_name <- colnames(data_cols)[2*i - 1]
      
      data.frame(
        Condition = meta[[1]],
        Gene = gene_name,
        E = e_col,
        Ct = ct_col,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  out[[1]] <- factor(out[[1]])
  out
}



.long_to_wide <- function(df) {
  
  if (ncol(df) < 4) {
    stop("Data frame must contain at least 4 columns: Condition, Gene, E, Ct")
  }
  
  # standardize first 4 columns internally
  tmp <- data.frame(
    Condition = df[[1]],
    Gene      = df[[2]],
    E         = df[[3]],
    Ct        = df[[4]],
    stringsAsFactors = FALSE
  )
  
  # replicate number within Condition C Gene
  tmp$Rep <- ave(
    seq_len(nrow(tmp)),
    tmp$Condition,
    tmp$Gene,
    FUN = seq_along
  )
  
  # reshape to wide
  wide <- reshape(
    tmp,
    idvar   = c("Condition", "Rep"),
    timevar = "Gene",
    direction = "wide"
  )
  wide[[1]] <- factor(wide[[1]])
  .cleanup(wide)
}



.geom_pub_cols <- function(col_width = 0.8,
                           err_width = 0.15,
                           fill_colors = NULL,
                           dodge_width = 0.8,
                           alpha = 1,
                           ...) {
  
  pos <- position_dodge(width = dodge_width)
  
  layers <- list(
    geom_col(width = col_width, position = pos, alpha = alpha, ...),
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = err_width, position = pos)
  )
  
  if (!is.null(fill_colors)) {
    layers <- c(layers, scale_fill_manual(values = fill_colors))
  }
  
  layers
}


.theme_pub <- function(base_size = 12,
                       base_family = "sans",
                       legend_position = "right",
                       ...) {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      legend.position = legend_position,
      ...
    )
}

# convert_to_character function
.convert_to_character <- function(numbers) {
  characters <- character(length(numbers))  # Initialize a character vector to store the results
  
  for (i in seq_along(numbers)) {
    if (numbers[i] < 0.001) {
      characters[i] <- "***"
    } else if (numbers[i] < 0.01) {
      characters[i] <- "**"
    } else if (numbers[i] < 0.05) {
      characters[i] <- "*"
    } else if (numbers[i] < 0.1) {
      characters[i] <- "."
    } else {
      characters[i] <- " "
    }
  }
  return(characters)
}



.cleanup <- function(x, numOfFactors, block) {
  
  # define factor columns 
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(x[seq_len(numOfFactors)], factor)
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(x[seq_len(numOfFactors + 1)], factor)
  }

  # pairwise cleanup (from last to first) 
  n <- ncol(x)
  
  for (i in seq(from = n, to = 2, by = -2)) {
    
    colA <- x[[i - 1]]
    colB <- x[[i]]
    
    # skip if either column is a factor
    if (is.factor(colA) || is.factor(colB)) next
    
    badA <- suppressWarnings(is.na(as.numeric(colA))) | colA == 0
    badB <- suppressWarnings(is.na(as.numeric(colB))) | colB == 0
    
    colA[badB] <- NA
    colB[badA] <- NA
    
    x[[i - 1]] <- colA
    x[[i]]     <- colB
  }

  
  x[] <- lapply(x, function(col) {
    if (is.factor(col)) return(col)
    if (is.character(col)) {
      col[col %in% c("Undetermined", "undetermined")] <- NA
      suppressWarnings(col <- as.numeric(col))
    }
    if (is.numeric(col)) {
      col[col == 0] <- NA
    }
    col
  })
  x
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

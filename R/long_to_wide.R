#' @title 
#' Converts a 4-column qPCR long data format to wide format
#' @description
#' Converts a 4-column (Condition, gene, Efficiency, Ct) qPCR long data format to wide format
#' @details
#' Converts a 4-column (Condition, gene, Efficiency, Ct) qPCR long data format to wide format
#' @author
#' Ghader Mirzaghaderi
#' @export
#' @param x 
#' a 4-column (Condition, gene, Efficiency, Ct) qPCR long data
#' @return
#'  A wide qPCR data frame
#'  
#' @examples 
#' 
#' df <- read.table(header = TRUE, text = "
#' Condition	Gene	E	Ct
#' control	C2H2-26	1.8	31.26
#' control	C2H2-26	1.8	31.01
#' control	C2H2-26	1.8	30.97
#' treatment	C2H2-26	1.8	32.65
#' treatment	C2H2-26	1.8	32.03
#' treatment	C2H2-26	1.8	32.4
#' control	C2H2-01	1.75	31.06
#' control	C2H2-01	1.75	30.41
#' control	C2H2-01	1.75	30.97
#' treatment	C2H2-01	1.75	28.85
#' treatment	C2H2-01	1.75	28.93
#' treatment	C2H2-01	1.75	28.9
#' control	C2H2-12	2	28.5
#' control	C2H2-12	2	28.4
#' control	C2H2-12	2	28.8
#' treatment	C2H2-12	2	27.9
#' treatment	C2H2-12	2	28
#' treatment	C2H2-12	2	27.9
#' control	ref	1.9	28.87
#' control	ref	1.9	28.42
#' control	ref	1.9	28.53
#' treatment	ref	1.9	28.31
#' treatment	ref	1.9	29.14
#' treatment	ref	1.9	28.63")
#' 
#' long_to_wide(df)


long_to_wide <- function(x) {
  
  numOfFactors <- 1
  block <- NULL
  df <- x
  
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
  
  # replicate number within Condition × Gene
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
  x <- wide
  
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
    
    # identify invalid cells
    badA <- suppressWarnings(is.na(as.numeric(colA))) | colA == 0
    badB <- suppressWarnings(is.na(as.numeric(colB))) | colB == 0
    
    # cross-invalidate within the pair
    colA[badB] <- NA
    colB[badA] <- NA
    
    x[[i - 1]] <- colA
    x[[i]]     <- colB
  }
  
  # column-wise cleanup (original behavior)
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

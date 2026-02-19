#' @title Bar plot of single gene expression (Delta Delta Ct method)
#' 
#' @description
#' Creates a bar plot of relative gene expression (fold change) values
#' from single gene analysis
#' 
#' @export
#' 
#' @param res An object created by the \code{ANOVA_DDCt()} function
#' @param col_width Numeric. Width of bars (default \code{0.8})
#' @param err_width Numeric. Width of error bars (default \code{0.15})
#' @param color Optional color for the bar outline
#' @param alpha Numeric. Transparency of bars (default \code{1})
#' @param base_size Numeric. Base font size for theme (default \code{12})
#' @param d Distance between horizontal significance lines
#' @param ... Additional ggplot2 layer arguments
#' 
#' @import ggplot2
#' @import ggsignif
#' 
#' @return
#' A bar plot of Delta Delta Ct showing pairwise significance 
#' 
#' @examples
#' res <- ANOVA_DDCt(
#' data_2factor,
#' numOfFactors = 2,
#' mainFactor.column = 2,
#' numberOfrefGenes = 1,
#' block = NULL,
#' analyseAllTarget = TRUE) # If you have multi-target gene data, specify a single target gene.
#' 
#' plotSingleGene(res, fill = "cyan4", color = "black", base_size = 12)

plotSingleGene <- function(res,
                           col_width = 0.8,
                           err_width = 0.15,
                           color = "black",
                           alpha = 1,
                           base_size = 12,
                           d = 0.4,
                           ...){
  
  if (length(unique(res$relativeExpression$gene)) == 1){
    gene_name <- unique(res$relativeExpression$gene)
  } else {
    stop("Only single gene analysis is allowed.")
  }
  
  
  lm <- res$perGene[[gene_name]]$lm
  df2 <- Means_DDCt(lm, specs = names(res$perGene[[gene_name]]$Final_data[1]))
  df1 <- res$relativeExpression
  df1$contrast <- sub(" .*", "", df1$contrast)
  df1$contrast <- factor(df1$contrast, levels = unique(df1$contrast))
  
  
  df2_filt <- df2[df2$sig != " ", ]
  df2_filt$contrast <- gsub(" vs ", ' ', df2_filt$contrast)
  df2_filt$contrast <- strsplit(df2_filt$contrast, " ")
  
  
  p <- ggplot(df1, aes(x = contrast, y = RE)) +
    geom_col(width = col_width, alpha = alpha, color = color, ...) +
    geom_errorbar(aes(ymin = Lower.se.RE, ymax = Upper.se.RE), width = err_width) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  if (nrow(df2_filt) > 0) {
    p <- p + geom_signif(comparisons = as.list(df2_filt$contrast),
                         y_position = seq(max(df1$Upper.se.RE),(max(df1$Upper.se.RE) + ((length(df1$contrast)*length(df1$contrast)))/2), by = d),
                         annotations = df2_filt$sig)
  }
  
  p + .theme_pub(base_size = base_size) +
    xlab(NULL) +
    theme(axis.text.x = element_text(size = base_size, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = base_size,color = "black", angle = 0),
          panel.border = element_rect(color = "black"))   # keeps the text visible
}
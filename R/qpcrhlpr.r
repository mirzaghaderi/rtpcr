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
    } else {
      characters[i] <- "ns"
    }
  }
  
  return(characters)
}


.invOrder <- function(invg){
    collapsed <- paste(invg,sep="", collapse = "")
    u <- unique(strsplit(collapsed, "")[[1]])
    if(length(u) < 2){
      return( invg)
    }
    u <- u[order(u)]
    m <- matrix(nrow = NROW(invg), ncol=length(u))
    m[] <-F
    for(i in 1:length( invg)){
      s <- strsplit( invg[i],"")[[1]]
      index <- match(s, u)
      m[i, index] <- T
    }
    for(i in 1:(length(u) - 1)){
      firstColT <- match(T, m[, i])[1]
      firstT <- match(T, rowSums(m[, i:length(u)] > 0))[1]
      if(firstT < firstColT){
        colT <- match(T, m[firstT, i:length(u)])[1]
        colT <- colT + i - 1
        tmp <- m[, colT]
        m[, colT] <- m[,i]
        m[, i] <- tmp
      }
    }
    res <- vector(mode = "character", length = length("trt"))
    for(i in 1:length(invg)){
      l <- u[m[i, ]]
      res[i] <- paste(l, sep = "",collapse = "")
    }
    return(res)
  }
  
.addwDCt <- function(x, 
                     block = NULL, 
                     numberOfrefGenes = 1){
  
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-5)]
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-7)]
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x[1:(ncol(x)-2)], wDCt = (log10(x$Etarget)*x$Cttarget)-
                        ((log10(x$Eref)*x$Ctref) + (log10(x$Eref2)*x$Ctref2))/2)
    }
    
  } else {
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-6)]
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      factors <- colnames(x)[1:(ncol(x)-8)]
      colnames(x)[ncol(x)-7] <- "block"
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x[1:(ncol(x)-2)], wDCt = (log10(x$Etarget)*x$Cttarget)-
                        ((log10(x$Eref)*x$Ctref) + (log10(x$Eref2)*x$Ctref2))/2)
    }
  }
  return(list(x = x, factors = factors))
}

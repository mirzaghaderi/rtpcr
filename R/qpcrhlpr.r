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

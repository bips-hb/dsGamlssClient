#'
#' @title Split input string by '$'
#' @description This is an internal function based on the internal \code{\link[dsBaseClient]{extract}} function from \code{dsBaseClient} (version 6.3.0). It splits
#' the input by the '$' symbol and returns the single elements as a list.
#' @param input A vector or a list of strings.
#' @keywords internal
#' @return a list with the following elements
#' \describe{
#' \item{\code{holders}}{The strings before the '$' symbol.}
#' \item{\code{elements}}{The strings after the '$' symbol.}
#' }
#' @author DataSHIELD Development Team
#'
extract <- function(input) {
  input <- unlist(input)
  output1 <- c()
  output2 <- c()
  for (i in 1:length(input)) {
    inputterms <- unlist(strsplit(input[i], "\\$", perl = TRUE))
    if (length(inputterms) > 1) {
      obj1 <- strsplit(input[i], "\\$", perl = TRUE)[[1]][1]
      obj2 <- strsplit(input[i], "\\$", perl = TRUE)[[1]][2]
    } else {
      obj1 <- NA
      obj2 <- strsplit(input[i], "\\$", perl = TRUE)[[1]][1]
    }
    output1 <- append(output1, obj1)
    output2 <- append(output2, obj2)
  }
  output <- list("holders" = output1, "elements" = output2)
  return(output)
}

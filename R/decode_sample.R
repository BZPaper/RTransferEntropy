#' Function to decode a combined and coded sample of two time series.
#'
#' @param x coded time series
#' @param coded table of codes and associated values
#'
#' @return returns a matrix
#' @export
#'
#' @examples
#' @keywords internal
#'
decode_sample <- function(x, coded) {

  n <- length(x)
  numvalues <- length(coded)
  index <- 1:n
  mat <- matrix(NA, ncol = 2, nrow = n)
  colnames(mat) <- c("x", "y")

  for(i in 1:numvalues) {
    pos <- index[x == coded[i]]
    if (length(pos) != 0) {
      mat[pos, 1] <- substr(names(coded[i]), 1, 1)
      mat[pos, 2] <- substr(names(coded[i]), 2, 2)
    }
  }

  return(mat)
}

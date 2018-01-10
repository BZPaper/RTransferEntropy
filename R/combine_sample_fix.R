#' Function that combines the given time series.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#'
#' @return returns a list
#' @export
#'
#' @examples
#'
combine_sample <- function(x, y) {

  n <- length(x)
  csample <- cbind(x, y)
  csample <- apply(csample, 1, function(x) paste(x, collapse = ""))

  valuestab <- sort(unique(csample))
  names(valuestab) <- valuestab

  return(list(valuestab = valuestab, csample = csample))
}

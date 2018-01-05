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
  newsample <- cbind(x, y)
  newsample <- apply(newsample, 1, function(x) paste(x, collapse = ""))

  values <- sort(unique(newsample))
  numvalues <- length(values)

  pool <- LETTERS
  valuestab <- pool[1:numvalues]
  names(valuestab) <- values

  csample   <- numeric(n)
  csample[] <- NA

  if (length(pool) < numvalues) {
    print("Too many combinations.")
    stop
  } else {
    for (i in 1:numvalues){
      csample <- ifelse(csample == values[i], pool[i], csample)
    }
  }

  return(list(valuestab = valuestab, csample = csample))
}

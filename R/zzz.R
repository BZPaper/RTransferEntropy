#' @useDynLib RTransferEntropy
#' @importFrom Rcpp sourceCpp
NULL


star <- function(x) {
  ifelse(is.null(x) || is.na(x), "",
         ifelse(x < 0.001, "***",
                ifelse(x < 0.01, "**",
                       ifelse(x < 0.05, "*",
                              ifelse(x < 0.1, ".", "")))))
}

fupper <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

#' Extract the Coefficient Matrix from a TEResult
#'
#' @param x a TEResult
#'
#' @return a Matrix containing the coefficients
#' @export
#'
#' @examples
#' set.seed(1234567890)
#' n <- 1000
#' x <- rep(0, n + 1)
#' y <- rep(0, n + 1)
#'
#' for (i in seq(n)) {
#'   x[i + 1] <- 0.2 * x[i] + rnorm(1, 0, 2)
#'   y[i + 1] <- x[i] + rnorm(1, 0, 2)
#' }
#'
#' x <- x[-1]
#' y <- y[-1]
#'
#' te_result <- transfer_entropy(x, y)
#' coefs(te_result)
coefs <- function(x) {
  if (!is.TEResult(x)) stop("x must be a TEResult")
  return(x$coef)
}

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
#' coef(te_result)
coef <- function(x) {
  if (!is.TEResult(x)) stop("x must be a TEResult")
  return(x$coef)
}

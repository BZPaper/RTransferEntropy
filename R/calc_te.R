#' Calculates the Transfer Entropy for two time series
#'
#' @inheritParams transfer_entropy
#'
#' @return a single numerical value for the transfer entropy
#' @export
#'
#' @seealso \code{\link{calc_ete}} and\code{\link{transfer_entropy}}
#' @examples
#' # construct two time-series
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
#' # calculate the X->Y transfer entropy value
#' calc_te(x, y)
#'
#' # calculate the Y->X transfer entropy value
#' calc_te(y, x)
#'
#' \donttest{
#'   # Compare the results
#'   calc_te(x, y, seed = 123)
#'   calc_te(y, x, seed = 123)
#'   transfer_entropy(x, y, nboot = 0, seed = 123)
#' }
calc_te <- function(x, y, lx = 1, ly = 1, q = 0.1,
                    entropy = "Shannon",
                    shuffles = 100,
                    type = "quantiles",
                    quantiles = c(5, 95),
                    bins = NULL,
                    limits = NULL,
                    burn = 50,
                    seed = NULL,
                    na.rm = TRUE) {
  calc_te_ete("te", x, y,
    lx = lx, ly = ly, entropy = entropy, q = q,
    shuffles = shuffles, type = type, quantiles = quantiles,
    bins = bins, limits = limits, burn = burn, seed = seed,
    na.rm = na.rm
  )
}

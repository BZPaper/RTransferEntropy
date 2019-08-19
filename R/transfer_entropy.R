#' Function to estimate Shannon and Renyi transfer entropy between two time
#' series x and y.
#'
#' @param x a vector of numeric values, ordered by time.
#'         Also allowed are \code{\link[xts]{xts}}, \code{\link[zoo]{zoo}},
#'         or \code{\link[stats]{ts}} objects.
#' @param y a vector of numeric values, ordered by time.
#'         Also allowed are \code{\link[xts]{xts}}, \code{\link[zoo]{zoo}},
#'         or \code{\link[stats]{ts}} objects.
#' @param lx Markov order of x, i.e. the number of lagged values affecting the
#'           current value of x. Default is \code{lx = 1}.
#' @param ly Markov order of y, i.e. the number of lagged values affecting the
#'           current value of y. Default is \code{ly = 1}.
#' @param q a weighting parameter used to estimate Renyi transfer entropy,
#'          parameter is between 0 and 1. For \code{q = 1}, Renyi transfer
#'          entropy converges to Shannon transfer entropy. Default is
#'          \code{q = 0.1}.
#' @param entropy specifies the transfer entropy measure that is estimated,
#'                either 'Shannon' or 'Renyi'. The first character can be used
#'                to specify the type of transfer entropy as well. Default is
#'                \code{entropy = 'Shannon'}.
#' @param shuffles the number of shuffles used to calculate the effective
#'                 transfer entropy. Default is \code{shuffles = 100}.
#' @param type specifies the type of discretization applied to the observed time
#'             series:'quantiles', 'bins' or 'limits'. Default is
#'             \code{type = 'quantiles'}.
#' @param quantiles specifies the quantiles of the empirical distribution of the
#'                  respective time series used for discretization.
#'                  Default is \code{quantiles = c(5,95)}.
#' @param bins specifies the number of bins with equal width used for
#'             discretization. Default is \code{bins = NULL}.
#' @param limits specifies the limits on values used for discretization.
#'               Default is \code{limits = NULL}.
#' @param nboot the number of bootstrap replications for each direction of
#'              the estimated transfer entropy. Default is \code{nboot = 300}.
#' @param burn the number of observations that are dropped from the beginning of
#'             the bootstrapped Markov chain. Default is \code{burn = 50}.
#' @param quiet if FALSE (default), the function gives feedback.
#' @param seed a seed that seeds the PRNG (will internally just call set.seed),
#'             default is \code{seed = NULL}.
#' @param na.rm if missing values should be removed (will remove the values at
#'             the same point in the other series as well). Default is \code{TRUE}.
#'
#' @return an object of class transfer_entropy, containing the transfer entropy
#'         estimates in both directions, the effective transfer entropy
#'         estimates in both directions, standard errors and p-values based on
#'         bootstrap replications of the Markov chains under the null hypothesis
#'         of statistical independence, an indication of statistical
#'         significance, and quantiles of the bootstrap samples
#'         (if \code{nboot > 0}).
#' @export
#'
#' @seealso \code{\link{coef}}, \code{\link{print.transfer_entropy}}
#'
#' @examples
#' # construct two time-series
#' set.seed(1234567890)
#' n <- 500
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
#' # Calculate Shannon's Transfer Entropy
#' te_result <- transfer_entropy(x, y, nboot = 100)
#' te_result
#'
#' summary(te_result)
#'
#' \donttest{
#'   # Parallel Processing using the future-package
#'   library(future)
#'   plan(multiprocess)
#'
#'   te_result2 <- transfer_entropy(x, y, nboot = 100)
#'   te_result2
#'
#'   # revert back to sequential execution
#'   plan(sequential)
#'
#'   te_result2 <- transfer_entropy(x, y, nboot = 100)
#'   te_result2
#'
#'   # General set of quiet
#'   set_quiet(TRUE)
#'   a <- transfer_entropy(x, y, nboot = 0)
#'
#'   set_quiet(FALSE)
#'   a <- transfer_entropy(x, y, nboot = 0)
#' }
transfer_entropy <- function(x,
                             y,
                             lx = 1,
                             ly = 1,
                             q = 0.1,
                             entropy = "Shannon",
                             shuffles = 100,
                             type = "quantiles",
                             quantiles = c(5, 95),
                             bins = NULL,
                             limits = NULL,
                             nboot = 300,
                             burn = 50,
                             quiet = NULL,
                             seed = NULL,
                             na.rm = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  t0 <- Sys.time()
  # Check for unequal length of time series
  if (length(x) != length(y)) {
    stop("x and y must be of same length.")
  }
  if (is.null(quiet)) quiet <- as.logical(options("RTransferEntropy::quiet"))

  # check and convert x and y if of class zoo (also includes xts)
  if ("zoo" %in% class(x)) x <- as.numeric(sort(x))
  if ("zoo" %in% class(y)) y <- as.numeric(sort(y))

  # Check that type is specified correctly
  type <- tolower(type)
  if (!type %in% c("quantiles", "bins", "limits", "q", "b", "l")) {
    stop("type must be either 'quantiles', 'bins' or 'limits'.")
  }

  if (nchar(type) == 1) {
    if (type == "q") {
      type <- "quantiles"
    } else if (type == "b") {
      type <- "bins"
    } else {
      type <- "limits"
    }
  }

  # Check/Restrict number of classes and Markov order/lags
  if (length(quantiles) > 20 || length(bins) > 20 || length(limits) > 20) {
    stop(paste(
      "Number of classes should not exceed 20.",
      "Do not expect sensical results when using too many classes and/or lags."
    ))
  }

  if (lx > 20 || ly > 20) {
    stop(paste(
      "Markov order/number of lags should not exceed 20.",
      "Do not expect sensical results when using too many classes and/or lags."
    ))
  }

  if (lx != ly) {
    warning(paste(
      "Markov order/number of lags should be identical for both",
      "time series to facilitate interpretation of results.",
      "Consider setting lx = ly."
    ))
  }

  # Check that transfer entropy measure is specified correctly
  entropy <- tolower(entropy)
  # Allow for specifying the first character only
  if (nchar(entropy) == 1 && entropy %in% c("s", "r")) {
    entropy <- if (entropy == "s") "shannon" else "renyi"
  }

  if (!entropy %in% c("shannon", "renyi")) {
    stop("entropy must be either 'Shannon' or 'Renyi'.")
  }

  # Check that q is between 0 and 1
  if (entropy == "renyi") {
    if (q < 0) {
      stop("q must follow 0 < q < 1")
    } else if (q >= 1) {
      warning(paste(
        "As q-->1, Renyi transfer entropy converges to Shannon transfer",
        "entropy. Shannon transfer entropy is calculated."
      ))
      entropy <- "shannon"
    }
  }

  # Check quantiles
  if (type == "quantiles" && (min(quantiles) < 0 || max(quantiles) > 100)) {
    stop("Quantiles must be between 0 and 100")
  }

  if (type == "quantiles" && max(quantiles) <= 1) {
    warning(paste(
      "Expected quantiles between 0 and 100 but found between 0 and 1,",
      "multiplying by 100."
    ))
    quantiles <- quantiles * 100
  }

  # Check number of bootstrap replications
  if (nboot > 0 && nboot < 100) {
    warning(paste(
      "Number of bootstrap replications is below 100. Use a higher number of",
      "bootstrap replications, you are relying on asymptotic arguments here."
    ))
  }

  x <- check_dimension(x)
  y <- check_dimension(y)

  # Remove missing values
  mis_values <- is.na(x) | is.na(y)
  if (na.rm == TRUE) {
    x <- x[!mis_values]
    y <- y[!mis_values]
  } else {
    if (any(mis_values)) return(NA)
  }

  if (length(x) == 0) stop("x and y must have non-missing values.")

  if (!quiet) {
    cat(sprintf(
      "%s's entropy on %s core%s with %s shuffle%s.\n",
      fupper(entropy),
      future::nbrOfWorkers(), mult_s(future::nbrOfWorkers()),
      shuffles, mult_s(shuffles)
    ))
    cat(sprintf(
      "  x and y have length %s (%s NAs removed)\n",
      length(x), sum(mis_values)
    ))
  }

  # Call either te_shannon or te_renyi
  if (entropy == "shannon") {
    te <- te_shannon(
      x = x,
      lx = lx,
      y = y,
      ly = ly,
      shuffles = shuffles,
      type = type,
      quantiles = quantiles,
      bins = bins,
      limits = limits,
      nboot = nboot,
      burn = burn,
      quiet = quiet
    )
  } else {
    te <- te_renyi(
      x = x,
      lx = lx,
      y = y,
      ly = ly,
      q = q,
      shuffles = shuffles,
      type = type,
      quantiles = quantiles,
      bins = bins,
      limits = limits,
      nboot = nboot,
      burn = burn,
      quiet = quiet
    )
  }


  # Inference (standard errors, p-values)
  if (nboot > 1) {
    seteyx <- sd(te$boot[1, ])
    setexy <- sd(te$boot[2, ])

    pval <- function(x, est) length(x[x > est]) / length(x)
    psteyx <- pval(te$boot[1, ], te$teyx)
    pstexy <- pval(te$boot[2, ], te$texy)
  } else {
    seteyx <- NA
    setexy <- NA

    psteyx <- NA
    pstexy <- NA
  }

  coef <- matrix(
    c(
      te$texy, te$stexy, setexy, pstexy,
      te$teyx, te$steyx, seteyx, psteyx
    ),
    nrow = 2, byrow = T,
    dimnames = list(
      c("X->Y", "Y->X"),
      c("te", "ete", "se", "p-value")
    )
  )
  if (entropy == "shannon") {
    coef[, "te"] <- pmax(0, coef[, "te"])
    coef[, "ete"] <- pmax(0, coef[, "ete"])
  }
  q <- ifelse(entropy == "renyi", q, NA)

  # Output
  res <- list(
    entropy = entropy,
    obs = list(x = x, y = y),
    coef = coef,
    nobs = length(x),
    q = q,
    boot = te$boot
  )

  if (!quiet) {
    t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat("Done - Total time", round(t, 2), "seconds\n")
  }

  class(res) <- "transfer_entropy"
  return(res)
}

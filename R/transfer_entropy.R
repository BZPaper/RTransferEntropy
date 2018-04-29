#' Calculates Shannon and Renyi transfer entropy between two time series.
#'
#' @param x vector of values, ordered by time
#' @param y vector of values, ordered by time
#' @param lx Markov order of x, i.e. number of lagged values affecting the
#'           current value; default is 1
#' @param ly Markov order of y, i.e. number of lagged values affecting the
#'           current value; default is 1
#' @param q weighting parameter in Renyi transfer entropy between 0 and 1;
#'          at \code{q = 1}, Renyi transfer entropy converges to Shannon
#'          transfer entropy; default is 0.1
#' @param entropy transfer entropy measure that is calculated, either 'Shannon'
#'                or 'Renyi'; first character can be used as well;
#'                default is Shannon
#' @param nreps number of replications for each shuffle; default is 2
#' @param shuffles number of shuffles; default is 50
#' @param cl numeric value (default is number of cores - 1),
#'           or a cluster as created by \code{\link[parallel]{makeCluster}}
#'           that can be used by \code{\link[pbapply]{pbapply}}
#' @param type 'quantiles', 'bins' or 'limits' to discretize the data; default
#'             is 'quantiles'
#' @param quantiles quantiles of empirical distribution used for discretization
#' @param bins number of bins with equal width used for discretization
#' @param limits user determined limits on values used for discretization
#' @param nboot number of bootstrap replications; default is 300
#' @param burn number of observations that are dropped from the beginning of
#'             the bootstrapped Markov chain; default is 50
#' @param quiet if FALSE (default), the function gives feedback
#' @param seed a seed that seeds the PRNG (will internally just call set.seed),
#'             default is NULL
#'
#' @return an object of class TEResult, containing the entropy measure, the
#'         effective transfer entropy measure, standard errors, p-values,
#'         indication of statistical significance, quantiles of bootstrap
#'         sample (if nboot > 0)
#' @export
#'
#' @seealso \code{\link{coef}}, \code{\link{print.TEResult}}
#'
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
#' te_result <- transfer_entropy(x, y)
#' te_result
#'
#' te_result <- transfer_entropy(x, y, nboot = 0)
#' te_result
#'
#' is.TEResult(te_result)
transfer_entropy <- function(x,
                             y,
                             lx = 1,
                             ly = 1,
                             q = 0.1,
                             entropy = "Shannon",
                             nreps = 2,
                             shuffles = 50,
                             cl = parallel::detectCores() - 1,
                             type = "quantiles",
                             quantiles = c(5, 95),
                             bins = NULL,
                             limits = NULL,
                             nboot = 300,
                             burn = 50,
                             quiet = FALSE,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  t0 <- Sys.time()
  # Check for unequal length of time series
  if (length(x) != length(y)) {
    stop("x and y must be of same length.")
  }

  # Check that type is specified correctly
  type <- tolower(type)
  if (!type %in% c("quantiles", "bins", "limits", "q", "b", "l"))
    stop("type must be either 'quantiles', 'bins' or 'limits'.")

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
  if (length(quantiles) > 20 || length(bins) > 20 || length(limits) > 20)
    stop("Number of classes should not exceed 20. Do not expect sensical results when using too many classes and/or lags.")

  if (lx > 20 || ly > 20)
    stop("Markov order/number of lags should not exceed 20. Do not expect sensical results when using too many classes and/or lags.")

  if (lx != ly)
    warning("Markov order/number of lags should be identical for both time series to facilitate interpretation of results. Consider setting lx = ly.")

  # Check that transfer entropy measure is specified correctly
  entropy <- tolower(entropy)
  # Allow for specifying the first character only
  if (nchar(entropy) == 1 && entropy %in% c("s", "r")) {
    entropy <- if (entropy == "s") "shannon" else "renyi"
  }

  if (!entropy %in% c("shannon", "renyi"))
    stop("entropy must be either 'Shannon' or 'Renyi'.")

  # Check that q is between 0 and 1
  if (entropy == "renyi") {
    if (q < 0) {
      stop("q must follow 0 < q < 1")
    } else if (q >= 1) {
      warning("As q-->1, Renyi transfer entropy converges to Shannon transfer entropy. Shannon transfer entropy is calculated.")
      entropy <- "shannon"
    }
  }

  # Check quantiles
  if (type == "quantiles" && (min(quantiles) < 0 || max(quantiles) > 100))
    stop("Quantiles must be between 0 and 100")

  if (type == "quantiles" && max(quantiles) <= 1) {
    warning("Expected quantiles between 0 and 100 but found between 0 and 1, multiplying by 100.")
    quantiles <- quantiles * 100
  }

  # Check number of bootstrap replications
  if (nboot < 100)
    warning("Number of bootstrap replications is below 100. Use a higher number of bootstrap replications, you are relying on asymptotic arguments here.")

  if (!quiet) cat(sprintf("Calculating %s's entropy ", fupper(entropy)))

  # Set-up the parallelization of computations
  if (is.numeric(cl)) {
    if (cl == 1) {
      cl <- NULL
      if (!quiet) cat("sequentially ")
    } else {
      cl <- min(cl, parallel::detectCores())
      if (!quiet) cat(sprintf("on %s cores ", cl))
      cl <- parallel::makeCluster(cl)
      on.exit(parallel::stopCluster(cl), add = T)
    }
  } else if ("cluster" %in% class(cl)) {
    if (!quiet) cat(sprintf("on %s cores ", length(cl)))
  } else  {
    stop("cl must be either a cluster (i.e., parallel::makeCluster()), or a numeric value")
  }

  # Only if a cluster is specified, set the pboptions
  # Otherwise, interference with potential pboptions set by the user
  if (!is.null(cl)) {
    pbapply::pboptions(type = "timer")
    if (quiet) pbapply::pboptions(type = "none")
  }

  # Remove missing values
  mis_values <- is.na(x) | is.na(y)
  x <- x[!mis_values]
  y <- y[!mis_values]

  if (!quiet) cat(sprintf("with %s shuffle(s) and %s bootstrap(s)\nThe timeseries have length %s (%s NAs removed)\n",
                          shuffles, nboot, length(x), sum(mis_values)))

  if (length(x) == 0) return(NA)

  # Call either te_shannon or te_renyi
  if (entropy == "shannon") {
    te <- te_shannon(x = x,
                     lx = lx,
                     y = y,
                     ly = ly,
                     nreps = nreps,
                     shuffles = shuffles,
                     cl = cl,
                     type = type,
                     quantiles = quantiles,
                     bins = bins,
                     limits = limits,
                     nboot = nboot,
                     burn = burn,
                     quiet = quiet)
  } else {
    te <- te_renyi(x = x,
                   lx = lx,
                   y = y,
                   ly = ly,
                   q = q,
                   nreps = nreps,
                   shuffles = shuffles,
                   cl = cl,
                   type = type,
                   quantiles = quantiles,
                   bins = bins,
                   limits = limits,
                   nboot = nboot,
                   burn = burn,
                   quiet = quiet)
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
    c(te$texy, max(0, te$stexy), setexy, pstexy,
      te$teyx, max(0, te$steyx), seteyx, psteyx),
    nrow = 2, byrow = T,
    dimnames = list(c("X->Y", "Y->X"),
                    c("te", "ete", "se", "p-value"))
  )

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

  class(res) <- append(class(res), "TEResult")
  return(res)
}

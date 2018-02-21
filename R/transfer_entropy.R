#' Wrapper for the implementation of Shannon and Renyi transfer entropy.
#'
#' @param x a vector of values
#' @param y a vector of values
#' @param lx Markov order of x, defaults to 1
#' @param ly Markov order of y, defaults to 1
#' @param q weighting parameter in Renyi transfer entropy, defaults to 0.5
#' @param entropy the type of entropy calculation to use, either 'shannon'
#'   or 'renyi', first character can be used as well, defaults to shannon.
#' @param constx constant value substracted from transfer entropy measure x,
#'   default to NULL (no const)
#' @param consty constant value substracted from transfer entropy measure y,
#'   default to NULL (no const)
#' @param nreps number of replications for each shuffle
#' @param shuffles number of shuffles
#' @param cl a numeric value (defaults to number of cores - 1),
#'    or a cluster as created by \code{\link[parallel]{makeCluster}}
#'    that can be used by \code{\link[pbapply]{pbapply}}
#' @param type bins, limits or quantiles of empirical distribution to discretize
#' the data
#' @param quantiles quantiles to use for discretization
#' @param bins the number of bins with equal width used for discretization
#' @param limits limits used for discretization
#' @param boots number of bootstrap samples
#' @param nboot number of bootstrap replications
#' @param burn number of observations that are dropped from the beginning of
#' the bootstrapped Markov chain
#' @param quiet if FALSE (default), the function gives feedback
#' @param seed a seed that seeds the PRNG (will internally just call set.seed), defaults to NULL
#'
#' @return an object of class TEResult, coontaining the entropy measure, the
#'   effective transfer entropy measure, standard errores, p-values, etc.
#' @export
#'
#' @examples
#' # construct two time-series
#' set.seed(1234567890)
#' n <- 100000
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
#' te_result <- transfer_entropy(x, lx = 1, y, ly = 1)
#' te_result
#'
#' summary(te_result)
#'
#' is.TEResult(te_result)
transfer_entropy <- function(x,
                             y,
                             lx = 1,
                             ly = 1,
                             q = 0.5,
                             entropy = "Shannon",
                             constx = NULL,
                             consty = NULL,
                             nreps = 2,
                             shuffles = 6,
                             cl = parallel::detectCores() - 1,
                             type = "quantiles",
                             quantiles = c(5, 95),
                             bins = NULL,
                             limits = NULL,
                             nboot = 10,
                             burn = 50,
                             quiet = FALSE,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  t0 <- Sys.time()
  # Check for unequal length of time series and treat missing values
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }

  # check that type is specified correctly
  type <- tolower(type)
  if (!type %in% c("quantiles", "bins", "limits", "q", "b", "l"))
    stop("type must be either 'quantiles', 'bins', or 'limits'")

  if (nchar(type) == 1) {
    if (type == "q") {
      type <- "quantiles"
    } else if (type == "b") {
      type <- "bins"
    } else {
      type <- "limits"
    }
  }

  # check that entropy is specified correctly
  entropy <- tolower(entropy)
  # allow to specify the first character only as well
  if (nchar(entropy) == 1 && entropy %in% c("s", "r")) {
    entropy <- if (entropy == "s") "shannon" else "renyi"
  }

  if (!entropy %in% c("shannon", "renyi"))
    stop("entropy must be either 'shannon' or 'renyi'.")

  if (!quiet) cat(sprintf("Calculating %s's entropy ", fupper(entropy)))

  # set-up the parallel stuff
  pbapply::pboptions(type = "timer")
  if (quiet) pbapply::pboptions(type = "none")

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
  } else {
    if (!quiet) cat(sprintf("on %s cores ", length(cl)))
  }

  # remove missing-values
  mis_values <- is.na(x) | is.na(y)
  x <- x[!mis_values]
  y <- y[!mis_values]

  if (!quiet) cat(sprintf("with %s shuffle(s) and %s bootstrap(s)\nThe timeseries have length %s (%s NAs removed)\n",
                          shuffles, nboot, length(x), sum(mis_values)))

  # call either te_shannon or te_renyi
  if (entropy == "shannon") {
    te <- te_shannon(x = x,
                     lx = lx,
                     y = y,
                     ly = ly,
                     constx = constx,
                     consty = consty,
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
                   constx = constx,
                   consty = consty,
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
    psteyx <- pval(te$boot[1, ], te$steyx)
    pstexy <- pval(te$boot[2, ], te$stexy)
  } else {
    seteyx <- NA
    setexy <- NA

    psteyx <- NA
    pstexy <- NA
  }

  coef <- matrix(
    c(te$texy, te$stexy, setexy, pstexy,
      te$teyx, te$steyx, seteyx, psteyx),
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

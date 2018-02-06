#' Wrapper for the implementation of Shannon and Renyi transfer entropy.
#'
#' @param x a vector of values
#' @param lx Markov order of x, defaults to 1
#' @param y a vector of values
#' @param ly Markov order of y, defaults to 1
#' @param q weighting parameter in Renyi transfer entropy, defaults to 0.5
#' @param entropy the type of entropy calculation to use, either 'shannon'
#'   or 'renyi', first character can be used as well, defaults to shannon.
#' @param shuffle if TRUE (default), shuffled transfer entropy is calculated
#' @param const if TRUE (defaults to FALSE), then shuffle is constant for all bootstraps
#' @param constx constant value substracted from transfer entropy measure x
#' @param consty constant value substracted from transfer entropy measure y
#' @param nreps number of replications for each shuffle
#' @param shuffles number of shuffles
#' @param ncores number of cores in parallel computation
#' @param type bins, limits or quantiles of empirical distribution to discretize
#' the data
#' @param quantiles quantiles to use for discretization
#' @param bins the number of bins with equal width used for discretization
#' @param limits limits used for discretization
#' @param boots number of bootstrap samples
#' @param nboot number of bootstrap replications
#' @param burn number of observations that are dropped from the beginning of
#' the bootstrapped Markov chain
#'
#' @return returns a list containing the respective transfer entropy measure,
#' the effective transfer entropy measure, standard errors, p-values and
#' number of observations
#' @export
#'
#' @examples
#'
transfer_entropy <- function(x,
                             lx = 1,
                             y,
                             ly = 1,
                             q = 0.5,
                             entropy = "Shannon",
                             shuffle = TRUE,
                             const = FALSE,
                             constx = 0,
                             consty = 0,
                             nreps = 2,
                             shuffles = 6,
                             ncores = parallel::detectCores() - 1,
                             type = "quantiles",
                             quantiles = c(5, 95),
                             bins = NULL,
                             limits = NULL,
                             nboot = 1000,
                             burn = 50) {

  # Check for unequal length of time series and treat missing values
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }

  tsmat <- na.omit(matrix(cbind(x, y), ncol = 2))
  obs <- dim(tsmat)[1]

  # Calculate transfer entropy
  entropy <- tolower(entropy)

  # allow to specify the first character only as well
  if (nchar(entropy) == 1 && entropy %in% c("s", "r")) {
    entropy <- if (entropy == "s") "shannon" else "renyi"
  }

  if (!entropy %in% c("shannon", "renyi"))
    stop("entropy must be either 'shannon' or 'renyi'.")

  # assign the respective function to te_function
  te_function <- if (entropy == "shannon") te_shannon else te_renyi

  # call either te_shannon or te_renyi as the te_function
  te <- te_function(x = tsmat[, 1],
                    lx = lx,
                    y = tsmat[, 2],
                    ly = ly,
                    shuffle = shuffle,
                    const = const,
                    constx = constx,
                    consty = consty,
                    nreps = nreps,
                    shuffles = shuffles,
                    ncores = ncores,
                    type = type,
                    quantiles = quantiles,
                    bins = bins,
                    limits = limits,
                    nboot = nboot,
                    burn = burn)

  if (is.null(dim(te$bootstrap_H0))) {
    # Inference (standard errors, p-values)
    setex <- 0
    setey <- 0

    pstex <- 0
    pstey <- 0
  } else {
    # Inference (standard errors, p-values)
    setex <- sd(te$bootstrap_H0[1, ])
    setey <- sd(te$bootstrap_H0[2, ])

    pval <- function(x, est) length(x[x > est]) / length(x)
    pstex <- pval(te$bootstrap_H0[1, ], te$stex)
    pstey <- pval(te$bootstrap_H0[2, ], te$stey)
  }

  # Output
  res <- list(TE_YX = te$tex,
              TE_XY = te$tey,
              ETE_YX = te$stex,
              ETE_XY = te$stey,
              SETE_YX = setex,
              SETE_XY = setey,
              PETE_YX = pstex,
              PETE_XY = pstey,
              tsobs = obs)

  return(res)
}

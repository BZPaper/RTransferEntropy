#' Wrapper for the implementation of Shannon and Renyi transfer entropy.
#' --------------------------PRELIMINARY-------------------------------
#'
#' @param x
#' @param lx
#' @param y
#' @param ly
#' @param q
#' @param entropy
#' @param shuffle
#' @param const
#' @param constx
#' @param consty
#' @param nreps
#' @param shuffles
#' @param ncores
#' @param type
#' @param quantiles
#' @param bins
#' @param limits
#' @param nboot
#' @param burn
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
  te_function <- if (entropy == "shannon") te_shannnon else te_renyi

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

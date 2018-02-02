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
TE <- function(x,
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
  m <- length(x)
  n <- length(y)

  if (m != n) {
    stop(cat(paste("Warning: Time series of unequal length",
                   "\n", "Execution halted")))
  }

  tsmat <- na.omit(matrix(cbind(x, y), ncol = 2))
  obs <- dim(tsmat)[1]

  # Calculate transfer entropy
  if (entropy == "Shannon") {
    te <- Shannon_transfer_entropy(x = tsmat[, 1],
                                   lx,
                                   y = tsmat[, 2],
                                   ly,
                                   shuffle,
                                   const,
                                   constx,
                                   consty,
                                   nreps,
                                   shuffles,
                                   ncores,
                                   type,
                                   quantiles,
                                   bins,
                                   limits,
                                   nboot,
                                   burn)
  } else if (entropy == "Renyi") {
    te <- Renyi_transfer_entropy(x = tsmat[, 1],
                                 lx,
                                 y = tsmat[, 2],
                                 ly,
                                 q,
                                 shuffle,
                                 const,
                                 constx,
                                 consty,
                                 nreps,
                                 shuffles,
                                 ncores,
                                 type,
                                 quantiles,
                                 bins,
                                 limits,
                                 nboot,
                                 burn)
  } else {
    stop(cat(paste("Warning: Transfer entropy measure not correctly specified",
                   "\n", "Execution halted")))
  }

  if (nboot != 0) {
    # Inference (standard errors, p-values)
    setex <- sd(te$bootstrap_H0[1, ])
    setey <- sd(te$bootstrap_H0[2, ])

    pval <- function(x, est) {length(x[x > est])/length(x)}
    pstex <- pval(te$bootstrap_H0[1, ], te$stex)
    pstey <- pval(te$bootstrap_H0[2, ], te$stey)
  } else {
    # Inference (standard errors, p-values)
    setex <- 0
    setey <- 0

    pstex <- 0
    pstey <- 0
  }

  # Output
  return(list(TE_YX = te$tex,
              TE_XY = te$tey,
              ETE_YX = te$stex,
              ETE_XY = te$stey,
              SETE_YX = setex,
              SETE_XY = setey,
              PETE_YX = pstex,
              PETE_XY = pstey,
              tsobs = obs))
}

#' Function for bootstrapping Shannon transfer entropy under H0 of independence
#' between time series x and y.
#'
#' @param x a vector of coded values
#' @param lx Markov order of x
#' @param y a vector of coded values
#' @param ly Markov order of y
#' @param burn number of observations that are dropped from the beginning of
#' the bootstrapped Markov chain
#' @param shuffle if TRUE, shuffled transfer entropy is calculated
#' @param const if TRUE, then shuffle is constant for all bootstraps
#' @param constx constant value substracted from transfer entropy measure
#' @param consty constant value substracted from transfer entropy measure
#' @param nreps number of replications for each shuffle
#' @param shuffles number of shuffles
#' @param ncores number of cores in parallel computation
#'
#' @return returns a vector
#' @export
#'
#' @examples
#'
trans_boot_H0 <- function(x,
                          lx,
                          y,
                          ly,
                          burn = 50,
                          shuffle = TRUE,
                          const = FALSE,
                          constx = NULL,
                          consty = NULL,
                          nreps = 2,
                          shuffles = 6,
                          ncores = parallel::detectCores() - 1) {

  bootx <- Markov_boot_step(x, lx, burn)
  booty <- Markov_boot_step(y, ly, burn)

  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- transfer_entropy_internal(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- transfer_entropy_internal(booty, lx = ly, x, ly = lx)$transentropy - consty
    } else {
      constx <- shuffled_transfer_entropy(bootx,
                                          lx = lx,
                                          y,
                                          ly = ly,
                                          nreps,
                                          shuffles,
                                          diff = TRUE,
                                          ncores)
      consty <- shuffled_transfer_entropy(booty,
                                          lx = ly,
                                          x,
                                          ly = lx,
                                          nreps,
                                          shuffles,
                                          diff = TRUE,
                                          ncores)

      # Lead = x
      dtex <- transfer_entropy_internal(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- transfer_entropy_internal(booty, lx = ly, x, ly = lx)$transentropy - consty
    }
  } else {
    # Lead = x
    dtex <- transfer_entropy_internal(bootx, lx = lx, y, ly = ly)$transentropy
    # Lead = y
    dtey <- transfer_entropy_internal(booty, lx = ly, x, ly = lx)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

#' Function for bootstrapping R?nyi transfer entropy under H0 of independence
#' between time series x and y.
#'
#' @param x a vector of coded values
#' @param lx Markov order of x
#' @param y a vector of coded values
#' @param ly Markov order of y
#' @param q weighting parameter in R?nyi transfer entropy
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
trans_boot_H0_ren <- function(x,
                              lx,
                              y,
                              ly,
                              q,
                              burn = 50,
                              shuffle = TRUE,
                              const = FALSE,
                              constx = NULL,
                              consty = NULL,
                              nreps = 2,
                              shuffles = 6,
                              ncores = parallel::detectCores() - 1) {

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy
      - constx
      # Lead = y
      dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy
      - consty
    } else {
      constx <- shuffled_transfer_entropy_ren(bootx,
                                              lx = lx,
                                              y,
                                              ly = ly,
                                              q,
                                              nreps,
                                              shuffles,
                                              diff = TRUE,
                                              ncores)
      consty <- shuffled_transfer_entropy_ren(booty,
                                              lx = ly,
                                              x,
                                              ly = lx,
                                              q,
                                              nreps,
                                              shuffles,
                                              diff = TRUE,
                                              ncores)

      # Lead = x
      dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy
      - constx
      # Lead = y
      dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy
      - consty
    }
  } else {
    # Lead = x
    dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy
    # Lead = y
    dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

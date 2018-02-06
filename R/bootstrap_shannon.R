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
#' @keywords internal
#' @export
#'
#' @examples
#'
bootstrap_shannon <- function(x,
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

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- calc_te_shannon(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- calc_te_shannon(booty, lx = ly, x, ly = lx)$transentropy - consty
    } else {
      constx <- shuffle_shannon(bootx,
                                lx = lx,
                                y,
                                ly = ly,
                                nreps,
                                shuffles,
                                diff = TRUE,
                                ncores)
      consty <- shuffle_shannon(booty,
                                lx = ly,
                                x,
                                ly = lx,
                                nreps,
                                shuffles,
                                diff = TRUE,
                                ncores)

      # Lead = x
      dtex <- calc_te_shannon(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- calc_te_shannon(booty, lx = ly, x, ly = lx)$transentropy - consty
    }
  } else {
    # Lead = x
    dtex <- calc_te_shannon(bootx, lx = lx, y, ly = ly)$transentropy
    # Lead = y
    dtey <- calc_te_shannon(booty, lx = ly, x, ly = lx)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

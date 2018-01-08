#' Function for bootstrapping transfer entropy under H0.
#'
#' @param x a vector of coded values
#' @param valuestab table of codes and associated values
#' @param shuffle if TRUE, shuffled transfer entropy is calculated
#' @param lx x(k)
#' @param ly y(j)
#' @param const if TRUE, then shuffle is constant for all bootstraps
#' @param constx
#' @param consty
#' @param nreps
#' @param shuffles
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
                          shuffle = TRUE,
                          const = FALSE,
                          constx = NULL,
                          consty = NULL,
                          nreps = 2,
                          shuffles = 6,
                          ncores = parallel::detectCores() - 1) {

  bootx <- Markov_boot_step(x, lx)
  booty <- Markov_boot_step(y, ly)

  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- transfer_entropy(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- transfer_entropy(booty, lx = ly, x, ly = lx)$transentropy - consty
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
      dtex <- transfer_entropy(bootx, lx = lx, y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- transfer_entropy(booty, lx = ly, x, ly = lx)$transentropy - consty
    }
  } else {
    # Lead = x
    dtex <- transfer_entropy(bootx, lx = lx, y, ly = ly)$transentropy
    # Lead = y
    dtey <- transfer_entropy(booty, lx = ly, x, ly = lx)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

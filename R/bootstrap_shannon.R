#' Function for bootstrapping Shannon transfer entropy under H0 of independence
#' between time series x and y.
#'
#' @inheritParams transfer_entropy
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

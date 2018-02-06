#' Function for bootstrapping R?nyi transfer entropy under H0 of independence
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
bootstrap_renyi <- function(x,
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
      dtex <- calc_te_renyi(bootx, lx = lx, y, ly = ly, q)$transentropy
      - constx
      # Lead = y
      dtey <- calc_te_renyi(booty, lx = ly, x, ly = lx, q)$transentropy
      - consty
    } else {
      constx <- shuffle_renyi(bootx,
                              lx = lx,
                              y,
                              ly = ly,
                              q,
                              nreps,
                              shuffles,
                              diff = TRUE,
                              ncores)
      consty <- shuffle_renyi(booty,
                              lx = ly,
                              x,
                              ly = lx,
                              q,
                              nreps,
                              shuffles,
                              diff = TRUE,
                              ncores)

      # Lead = x
      dtex <- calc_te_renyi(bootx, lx = lx, y, ly = ly, q)$transentropy
      - constx
      # Lead = y
      dtey <- calc_te_renyi(booty, lx = ly, x, ly = lx, q)$transentropy
      - consty
    }
  } else {
    # Lead = x
    dtex <- calc_te_renyi(bootx, lx = lx, y, ly = ly, q)$transentropy
    # Lead = y
    dtey <- calc_te_renyi(booty, lx = ly, x, ly = lx, q)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

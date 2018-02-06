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
                              constx = NULL,
                              consty = NULL,
                              nreps = 2,
                              shuffles = 6,
                              ncores = parallel::detectCores() - 1) {

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  if (shuffle) {
    if (is.null(constx) || is.null(consty)) {
      # Lead = x
      dtex <- calc_te_shannon(x = bootx, lx = lx, y = y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- calc_te_shannon(x = booty, lx = ly, y = x, ly = lx)$transentropy - consty
    } else {
      constx <- shuffle_shannon(x = bootx,
                                lx = lx,
                                y = y,
                                ly = ly,
                                nreps = nreps,
                                shuffles = shuffles,
                                diff = TRUE,
                                ncores = ncores)

      consty <- shuffle_shannon(x = booty,
                                lx = ly,
                                y = x,
                                ly = lx,
                                nreps = nreps,
                                shuffles = shuffles,
                                diff = TRUE,
                                ncores = ncores)

      # Lead = x
      dtex <- calc_te_shannon(x = bootx, lx = lx, y = y, ly = ly)$transentropy - constx
      # Lead = y
      dtey <- calc_te_shannon(x = booty, lx = ly, y = x, ly = lx)$transentropy - consty
    }
  } else {
    # Lead = x
    dtex <- calc_te_shannon(x = bootx, lx = lx, y = y, ly = ly)$transentropy
    # Lead = y
    dtey <- calc_te_shannon(x = booty, lx = ly, y = x, ly = lx)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

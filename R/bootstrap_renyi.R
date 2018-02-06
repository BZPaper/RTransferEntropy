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
                            constx = 0,
                            consty = 0,
                            nreps = 2,
                            shuffles = 6,
                            ncores = parallel::detectCores() - 1) {

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  if (shuffle) {
    if (is.null(constx) || is.null(consty)) {
      # Lead = x
      x_te <- calc_te_renyi(x = bootx,
                            lx = lx, 
                            y = y, 
                            ly = ly, 
                            q = q)

      dtex <- x_te$transentropy - constx
      
      # Lead = y
      y_te <- calc_te_renyi(x = booty, 
                            lx = ly, 
                            y = x, 
                            ly = lx, 
                            q = q)

      dtey <- y_te$transentropy - consty
    } else {
      constx <- shuffle_renyi(x = bootx,
                              lx = lx,
                              y = y,
                              ly = ly,
                              q = q,
                              nreps = nreps,
                              shuffles = shuffles,
                              diff = TRUE,
                              ncores)
      
      consty <- shuffle_renyi(x = booty,
                              lx = ly,
                              y = x,
                              ly = lx,
                              q = q,
                              nreps = nreps,
                              shuffles = shuffles,
                              diff = TRUE,
                              ncores = ncores)

      # Lead = x
      x_te <- calc_te_renyi(x = bootx, 
                            lx = lx, 
                            y = y, 
                            ly = ly, 
                            q = q)

      dtex <- x_te$transentropy - constx
      
      # Lead = y
      y_te <- calc_te_renyi(x = booty, 
                            lx = ly, 
                            y = x, 
                            ly = lx, 
                            q = q)

      dtey <- y_te$transentropy - consty
    }
  } else {
    # Lead = x
    dtex <- calc_te_renyi(x = bootx, lx = lx, y = y, ly = ly, q = q)$transentropy
    # Lead = y
    dtey <- calc_te_renyi(x = booty, lx = ly, y = x, ly = lx, q = q)$transentropy
  }

  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")

  return(teboot)
}

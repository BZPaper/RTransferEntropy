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
                            constx = 0,
                            consty = 0,
                            nreps = 2,
                            shuffles = 6,
                            cl = NULL) {

  if (is.null(cl[[1]])) {
    opb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(opb), add = T)
  }

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  if (shuffles > 1) {
    if (is.null(constx) || is.null(consty)) {
      # Lead = x
      x_te <- calc_te_renyi(x = bootx,
                            lx = lx,
                            y = y,
                            ly = ly,
                            q = q)

      dteyx <- x_te$transentropy - constx

      # Lead = y
      y_te <- calc_te_renyi(x = booty,
                            lx = ly,
                            y = x,
                            ly = lx,
                            q = q)

      dtexy <- y_te$transentropy - consty
    } else {
      constx <- shuffle_renyi(x = bootx,
                              lx = lx,
                              y = y,
                              ly = ly,
                              q = q,
                              nreps = nreps,
                              shuffles = shuffles,
                              diff = TRUE,
                              cl = cl)

      consty <- shuffle_renyi(x = booty,
                              lx = ly,
                              y = x,
                              ly = lx,
                              q = q,
                              nreps = nreps,
                              shuffles = shuffles,
                              diff = TRUE,
                              cl = cl)

      # Lead = x
      x_te <- calc_te_renyi(x = bootx,
                            lx = lx,
                            y = y,
                            ly = ly,
                            q = q)

      dteyx <- x_te$transentropy - constx

      # Lead = y
      y_te <- calc_te_renyi(x = booty,
                            lx = ly,
                            y = x,
                            ly = lx,
                            q = q)

      dtexy <- y_te$transentropy - consty
    }
  } else {
    # Lead = x
    dteyx <- calc_te_renyi(x = bootx, lx = lx, y = y, ly = ly, q = q)$transentropy
    # Lead = y
    dtexy <- calc_te_renyi(x = booty, lx = ly, y = x, ly = lx, q = q)$transentropy
  }

  teboot <- c(dteyx, dtexy)
  names(teboot) <- c("dteyx", "dtexy")

  return(teboot)
}

#' Function to implement Shannon transfer entropy.
#'
#' @inheritParams transfer_entropy
#'
#' @return returns a list
#' @keywords internal
#' @export
#'
#' @examples
#'
te_shannon <- function(x,
                       lx,
                       y,
                       ly,
                       const = FALSE,
                       constx = 0,
                       consty = 0,
                       nreps = 2,
                       shuffles = 6,
                       cl = NULL,
                       type = "quantiles",
                       quantiles = c(5, 95),
                       bins = NULL,
                       limits = NULL,
                       nboot,
                       burn = 50,
                       quiet = FALSE) {

  # Code time series
  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # Lead = x
  if (!quiet) cat("Calculate the x->y transfer entropy\n")
  tex <- calc_te_shannon(x, lx = lx, y, ly = ly)$transentropy
  constx <- shuffle_shannon(x = x,
                            lx = lx,
                            y = y,
                            ly = ly,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            cl = cl)

  stex <- tex - constx

  # Lead = y
  if (!quiet) cat("Calculate the y->x transfer entropy\n")
  tey <- calc_te_shannon(y, lx = ly, x, ly = lx)$transentropy
  consty <- shuffle_shannon(x = y,
                            lx = ly,
                            y = x,
                            ly = lx,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            cl = cl)

  stey <- tey - consty

  # Bootstrap
  if (!quiet) cat("Bootstrap the transfer entropy\n")

  boot <- pbapply::pbreplicate(nboot,
                               bootstrap_shannon(x = x,
                                                 lx = lx,
                                                 y = y,
                                                 ly = ly,
                                                 burn = burn,
                                                 constx = constx,
                                                 consty = consty,
                                                 nreps = nreps,
                                                 shuffles = shuffles,
                                                 cl = NULL),
                               cl = cl)

  return(list(tex  = tex,
              tey  = tey,
              stex = stex,
              stey = stey,
              bootstrap_H0 = boot))
}

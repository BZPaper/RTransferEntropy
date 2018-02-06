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
                       shuffle = TRUE,
                       const = FALSE,
                       constx = 0,
                       consty = 0,
                       nreps = 2,
                       shuffles = 6,
                       ncores = parallel::detectCores() - 1,
                       type = "quantiles",
                       quantiles = c(5, 95),
                       bins = NULL,
                       limits = NULL,
                       nboot,
                       burn = 50) {

  # Code time series
  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # Calculate transfer entropy (without shuffling)
  # Lead = x
  tex <- calc_te_shannon(x, lx = lx, y, ly = ly)$transentropy
  # Lead = y
  tey <- calc_te_shannon(y, lx = ly, x, ly = lx)$transentropy

  # Calculate transfer entropy (with shuffling)
  constx <- shuffle_shannon(x = x,
                            lx = lx,
                            y = y,
                            ly = ly,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            ncores = ncores)

  consty <- shuffle_shannon(x = y,
                            lx = ly,
                            y = x,
                            ly = lx,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            ncores = ncores)

  # Lead = x
  stex <- tex - constx
  # Lead = y
  stey <- tey - consty

  # Bootstrap
  boot <- replicate(nboot,
                    bootstrap_shannon(x,
                                      lx = lx,
                                      y,
                                      ly = ly,
                                      burn,
                                      shuffle,
                                      const,
                                      constx,
                                      consty,
                                      nreps,
                                      shuffles,
                                      ncores))

  return(list(tex   = tex,
              tey   = tey,
              stex = stex,
              stey = stey,
              bootstrap_H0 = boot))
}

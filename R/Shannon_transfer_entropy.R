#' Function to implement Shannon transfer entropy.
#'
#' @param x
#' @param y
#' @param lx
#' @param ly
#' @param nreps
#' @param shuffles
#' @param const
#' @param shuffle
#' @param bins
#' @param quantiles
#' @param nboot
#'
#' @return returns a list
#' @export
#'
#' @examples
#'
Shannon_transfer_entopy <- function(x,
                                    lx,
                                    y,
                                    ly,
                                    shuffle = TRUE,
                                    const = FALSE,
                                    nreps = 2,
                                    shuffles = 6,
                                    ncores = parallel::detectCores() - 1,
                                    quantiles = c(5, 95),
                                    bins = NULL,
                                    limits = NULL,
                                    nboot) {

  # Code time series
  x <- code_sample(x, type = "quantiles", quantiles, bins, limits)
  y <- code_sample(y, type = "quantiles", quantiles, bins, limits)

  # Calculate transfer entropy (without shuffling)
  # Lead = x
  tex <- transfer_entropy(x, lx = lx, y, ly = ly)$transentropy
  # Lead = y
  tey <- transfer_entropy(y, lx = ly, x, ly = lx)$transentropy

  # Calculate transfer entropy (with shuffling)
  constx <- shuffled_transfer_entropy(x,
                                      lx = lx,
                                      y,
                                      ly = ly,
                                      nreps,
                                      shuffles,
                                      diff = FALSE,
                                      ncores)

  consty <- shuffled_transfer_entropy(y,
                                      lx = ly,
                                      x,
                                      ly = lx,
                                      nreps,
                                      shuffles,
                                      diff = FALSE,
                                      ncores)

  # Lead = x
  stex <- tex - constx
  # Lead = y
  stey <- tey - consty

  # Bootstrap
  boot1 <- replicate(nboot,
                     trans_boot_H0(x,
                                   lx = lx,
                                   y,
                                   ly = ly,
                                   shuffle,
                                   const = TRUE,
                                   constx = 0,
                                   consty = 0,
                                   nreps,
                                   shuffles,
                                   ncores))

  # Combine sample
  collapse <- combine_sample(x, y)
  collxy <- collapse$csample

  # Bootstrap
  boot2 <- replicate(nboot,
                     trans_boot_H1(collxy,
                                   lx = lx,
                                   ly = ly,
                                   collapse$valuestab,
                                   shuffle,
                                   const = TRUE,
                                   constx = mean(boot["dtex",]),
                                   consty = mean(boot["dtey",]),
                                   nreps,
                                   shuffles,
                                   ncores))

  return(list(tex   = tex,
              tey   = tey,
              S_tex = S_tex,
              S_tey = S_tey,
              bootstrap_H0 = boot1))
}

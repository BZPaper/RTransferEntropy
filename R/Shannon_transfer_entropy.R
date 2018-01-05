#' Function to implement Shannon transfer entropy.
#'
#' @param x
#' @param y
#' @param lx
#' @param ly
#' @param nrep
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
                                    y,
                                    lx,
                                    ly,
                                    nrep = 10,
                                    shuffles = 6,
                                    const = FALSE,
                                    shuffle = TRUE,
                                    bins,
                                    quantiles,
                                    nboot) {

  # Code time series
  x <- code_sample(x, type = "quantiles", quantiles)
  y <- code_sample(y, type = "quantiles", quantiles)

  # Calculate transfer entropy (withour shuffling)
  # Lead = x
  tex <- transfer_entropy(x, y, lx = lx, ly = ly)$transentropy
  # Lead = y
  tey <- transfer_entropy(y, x, lx = ly, ly = lx)$transentropy

  # Calculate transfer entropy (withour shuffling)
  constx <- shuffled_transfer_entropy(nrep,
                                      shuffles,
                                      diff = FALSE,
                                      x,
                                      lx = lx,
                                      y,
                                      ly = ly)

  consty <- shuffled_transfer_entropy(nrep,
                                      shuffles,
                                      diff = FALSE,
                                      y,
                                      lx = ly,
                                      x,
                                      ly = lx)

  # Lead = x
  stex <- tex - constx
  # Lead = y
  stey <- tey - consty

  # Bootstrap
  boot1 <- replicate(nboot,
                     trans_boot_H0(x,
                                   y,
                                   shuffle,
                                   lx = lx,
                                   ly = ly,
                                   const = TRUE,
                                   constx = 0,
                                   consty = 0,
                                   nrep,
                                   shuffles))

  # Combine sample
  collapse <- combine_sample(x, y)
  collxy <- collapse$csample

  # Calculate frequencies
  boot2 <- replicate(nboot,
                     trans_boot_H1(collxy,
                                   collapse$valuestab,
                                   shuffle,
                                   lx = lx,
                                   ly = ly,
                                   const = TRUE,
                                   constx = mean(boot["dtex",]),
                                   consty = mean(boot["dtey",]),
                                   nrep,
                                   shuffles))

  return(list(tex   = tex ,
              tey   = tey,
              S_tex = S_tex,
              S_tey = S_tey,
              bootstrap_H0 = boot1))
}

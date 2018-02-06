#' Function to implement Renyi transfer entropy.
#'
#' @inheritParams transfer_entropy
#'
#' @return returns a list
#' @keywords internal
#' @export
#'
#' @examples
#'
te_renyi <- function(x,
                     lx,
                     y,
                     ly,
                     q,
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

  # Calculate transfer entropy (withour shuffling)
  # Lead = x
  tex <- calc_te_renyi(x, lx = lx, y, ly = ly, q)$transentropy
  # Lead = y
  tey <- calc_te_renyi(y, lx = ly, x, ly = lx, q)$transentropy

  # Calculate transfer entropy (with shuffling)
  constx <- shuffle_renyi(x = x,
                          lx = lx,
                          y = y,
                          ly = ly,
                          q = q,
                          nreps = nreps,
                          shuffles = shuffles,
                          diff = FALSE,
                          ncores = ncores)

  consty <- shuffle_renyi(x = y,
                          lx = ly,
                          y = x,
                          ly = lx,
                          q = q,
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
                    bootstrap_renyi(x,
                                    lx,
                                    y,
                                    ly,
                                    q,
                                    burn,
                                    shuffle,
                                    const,
                                    constx,
                                    consty,
                                    nreps))


  return(list(tex   = tex,
              tey   = tey,
              stex = stex,
              stey = stey,
              bootstrap_H0 = boot))
}

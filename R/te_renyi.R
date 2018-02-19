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
                     nboot = 3,
                     burn = 50,
                     quiet = FALSE) {

  # Code time series
  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # Lead = x
  if (!quiet) cat("Calculate the x->y transfer entropy\n")
  tex <- calc_te_renyi(x, lx = lx, y, ly = ly, q)$transentropy
  constx <- shuffle_renyi(x = x,
                          lx = lx,
                          y = y,
                          ly = ly,
                          q = q,
                          nreps = nreps,
                          shuffles = shuffles,
                          diff = FALSE,
                          cl = cl)

  stex <- tex - constx

  # Lead = y
  if (!quiet) cat("Calculate the y->x transfer entropy\n")
  tey <- calc_te_renyi(y, lx = ly, x, ly = lx, q)$transentropy
  consty <- shuffle_renyi(x = y,
                          lx = ly,
                          y = x,
                          ly = lx,
                          q = q,
                          nreps = nreps,
                          shuffles = shuffles,
                          diff = FALSE,
                          cl = cl)

  stey <- tey - consty

  # Bootstrap
  if (!quiet) cat("Bootstrap the transfer entropy\n")

  seeds <- sample(.Machine$integer.max, nboot)

  boot <- pbapply::pbsapply(seeds, function(seed) {
    set.seed(seed)
    bootstrap_renyi(x = x,
                    lx = lx,
                    y = y,
                    ly = ly,
                    q = q,
                    burn = burn,
                    shuffles = shuffles,
                    constx = constx,
                    consty = consty,
                    nreps = nreps,
                    cl = NULL)
  }, cl = cl)

  return(list(tex   = tex,
              tey   = tey,
              stex = stex,
              stey = stey,
              bootstrap_H0 = boot))
}

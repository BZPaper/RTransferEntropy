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

  # Lead = y
  if (!quiet) cat("Calculate the X->Y transfer entropy\n")
  texy <- calc_te_shannon(y, lx = ly, x, ly = lx)$transentropy
  consty <- shuffle_shannon(x = y,
                            lx = ly,
                            y = x,
                            ly = lx,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            cl = cl)
  stexy <- texy - consty

  # Lead = x
  if (!quiet) cat("Calculate the Y->X transfer entropy\n")
  teyx <- calc_te_shannon(x, lx = lx, y, ly = ly)$transentropy
  constx <- shuffle_shannon(x = x,
                            lx = lx,
                            y = y,
                            ly = ly,
                            nreps = nreps,
                            shuffles = shuffles,
                            diff = FALSE,
                            cl = cl)
  steyx <- teyx - constx

  # Bootstrap
  if (!quiet) cat("Bootstrap the transfer entropies\n")
  if (nboot > 0) {
    seeds <- sample(.Machine$integer.max, nboot)
    boot <- pbapply::pbsapply(seeds, function(seed) {
      set.seed(seed)
      bootstrap_shannon(x = x,
                        lx = lx,
                        y = y,
                        ly = ly,
                        burn = burn,
                        constx = constx,
                        consty = consty,
                        nreps = nreps,
                        shuffles = shuffles,
                        cl = NULL)
    }, cl = cl)
  } else {
    boot <- NA
  }

  return(list(teyx  = teyx,
              texy  = texy,
              steyx = steyx,
              stexy = stexy,
              boot = boot))
}

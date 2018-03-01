# Function to implement Renyi transfer entropy.
# Same arguments as transfer_entropy.
# Used internally by transfer_entropy.
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

  # Lead = y
  if (!quiet) cat("Calculate the X->Y transfer entropy\n")
  texy <- calc_te_renyi(y, lx = ly, x, ly = lx, q)$transentropy
  consty <- shuffle_renyi(x = y,
                          lx = ly,
                          y = x,
                          ly = lx,
                          q = q,
                          nreps = nreps,
                          shuffles = shuffles,
                          diff = FALSE,
                          cl = cl)
  stexy <- texy - consty

  # Lead = x
  if (!quiet) cat("Calculate the Y->X transfer entropy\n")
  teyx <- calc_te_renyi(x, lx = lx, y, ly = ly, q)$transentropy
  constx <- shuffle_renyi(x = x,
                          lx = lx,
                          y = y,
                          ly = ly,
                          q = q,
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
  } else {
    boot <- NA
  }

  return(list(teyx   = teyx,
              texy   = texy,
              steyx = steyx,
              stexy = stexy,
              boot = boot))
}

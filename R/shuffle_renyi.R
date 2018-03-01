# Function to calculate the effective transfer entropy, as given by the
# difference between the Renyi transfer entropy calculated from a sample and
# the respective shuffled transfer entropy.
#
# Same arguments as transfer_entropy.
# Used internally by transfer_entropy.
# @return returns a numeric scalar
#
shuffle_renyi <- function(x,
                          lx,
                          y,
                          ly,
                          q,
                          nreps = 2,
                          shuffles = 6,
                          diff = TRUE,
                          cl = NULL) {

  seeds <- sample(.Machine$integer.max, shuffles)
  n <- length(x)

  shuffle <- pbapply::pblapply(seeds, function(seed) {
    set.seed(seed)
    res <- replicate(nreps,
                     calc_te_renyi(x = x,
                                   y = sample(y, n, replace = TRUE),
                                   lx = lx, ly = ly, q = q)$transentropy)
    return(res)
  }, cl = cl)

  ste <- mean(unlist(shuffle))

  if (diff) {
    te <- calc_te_renyi(x = x, y = y, lx = lx, ly = ly, q)$transentropy - ste
  } else {
    te <- ste
  }

  return(te)
}

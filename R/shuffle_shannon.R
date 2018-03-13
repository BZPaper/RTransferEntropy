# Function to calculate the effective transfer entropy, as given by the
# difference between the Shannon transfer entropy calculated from a sample and
# the respective shuffled transfer entropy. Used internally by transfer_entropy;
# same arguments.
#
shuffle_shannon <- function(x,
                            lx,
                            y,
                            ly,
                            nreps = 2,
                            shuffles = 6,
                            diff = TRUE,
                            cl = NULL) {

  seeds <- sample(.Machine$integer.max, shuffles)
  n <- length(x)

  shuffle <- pbapply::pblapply(seeds, function(seed) {
    set.seed(seed)
    res <- replicate(nreps,
                     calc_te_shannon(x = x,
                                     y = sample(y, n, replace = TRUE),
                                     lx = lx, ly = ly))
    return(res)
  }, cl = cl)

  ste <- mean(unlist(shuffle))

  if (diff) {
    te <- calc_te_shannon(x = x, y = y, lx = lx, ly = ly) - ste
  } else {
    te <- ste
  }

  return(te)
}

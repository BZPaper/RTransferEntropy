# Function to calculate the effective transfer entropy, as given by the
# difference between the Shannon transfer entropy calculated from a sample and
# the respective shuffled transfer entropy. Used internally by transfer_entropy;
# same arguments.
#
shuffle_shannon <- function(x, lx, y, ly, shuffles) {

  n <- length(x)

  if (shuffles > 200) {
    shuffle <- future.apply::future_replicate(
      shuffles,
      calc_te_shannon(x = x, y = sample(y, n, replace = TRUE), lx = lx, ly = ly)
    )
  } else {
    shuffle <- replicate(
      shuffles,
      calc_te_shannon(x = x, y = sample(y, n, replace = TRUE), lx = lx, ly = ly)
    )
  }

  ste <- mean(shuffle)

  return(ste)
}

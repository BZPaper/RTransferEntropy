# Function to calculate the effective transfer entropy, as given by the
# difference between the Renyi transfer entropy calculated from a sample and
# the respective shuffled transfer entropy. Used internally by transfer_entropy;
# same arguments.
#
shuffle_renyi <- function(x, lx, y, ly, q, shuffles) {
  n <- length(x)

  if (shuffles > 200) {
    shuffle <- future.apply::future_replicate(
      shuffles,
      calc_te_renyi(
        x = x, y = sample(y, n, replace = TRUE), lx = lx, ly = ly,
        q = q
      )
    )
  } else {
    shuffle <- replicate(
      shuffles,
      calc_te_renyi(
        x = x, y = sample(y, n, replace = TRUE), lx = lx, ly = ly,
        q = q
      )
    )
  }

  ste <- mean(shuffle)

  return(ste)
}

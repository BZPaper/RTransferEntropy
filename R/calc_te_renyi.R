# Function that calculates the transfer entropy between two time series x and
# y. The information flow from y to x is measured. Change x and y in function
# call to infer the dominant direction of the information flow. Calculated
# transfer entropy measure is Renyi transfer entropy. Used internally by
# transfer_entropy; same arguments.
#
calc_te_renyi <- function(x, lx, y, ly, q) {

  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, lx = lx, y, ly = ly)
  k1_j <- k1_j ^ q

  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)
  k1 <- k1 ^ q

  # x(k) and y(j)
  k_j <- cluster_gen(x, lx = lx, y, ly = ly, prog = FALSE)
  k_j <- k_j ^ q

  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)
  k <- k ^ q

  # Renyi transfer entropy
  #------------------------------

  numerator <- k1 / sum(k)
  denominator <- k1_j / sum(k_j)

  ren_entropy <- 1 / (1 - q) * log2(sum(numerator) / sum(denominator))

  return(ren_entropy)
}

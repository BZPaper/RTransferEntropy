# Function that calculates the transfer entropy between two time series x and
# y. The information flow from y to x is measured. Change x and y in function
# call to infer the dominant direction of the information flow. Calculated
# transfer entropy measure is Shannon transfer entropy. Used internally by
# transfer_entropy; same arguments.
#
calc_te_shannon <- function(x, lx, y, ly) {

  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, lx = lx, y, ly = ly)
  nck1_j <- length(k1_j)

  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)

  # x(k) and y(j)
  k_j <- cluster_gen(x, lx = lx, y, ly = ly, prog = FALSE)

  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)

  # Transfer entropy
  #------------------------------
  entropy <- numeric(nck1_j)
  for (i in 1:nck1_j) {
    names_ <- strsplit(names(k1_j[i]), " ")[[1]]

    p1 <- k1[paste0(names_[1:(lx + 1)], collapse = " ")]
    p2 <- k_j[paste0(names_[-(lx + 1)], collapse = " ")]
    p3 <- k[paste0(names_[1:lx], collapse = " ")]
    entropy[i] <- k1_j[i] * log2((k1_j[i] * p3) / (p2 * p1))
  }

  shan_entropy <- sum(entropy)

  return(shan_entropy)
}

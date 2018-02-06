#' Function that calculates the transfer entropy between two time series x and
#' y. The information flow from y to x is measured. Change x and y in function
#' call to infer the dominant direction of the information flow. Calculated
#' transfer entropy measure is Renyi transfer entropy.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#' @param lx Markov order of x
#' @param ly Markov order of y
#' @param q weighting parameter in R?nyi transfer entropy
#'
#' @keywords internal
#' @return returns a list
#' @export
#'
#' @examples
#'
calc_te_renyi <- function(x,
                          lx,
                          y,
                          ly,
                          q) {

  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, lx = lx, y, ly = ly)$frequency
  k1_j <- k1_j^q
  nck1_j <- length(k1_j)

  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)$frequency
  k1 <- k1^q
  nck1 <- length(k1)

  # x(k) and y(j)
  k_j <- cluster_gen(x, lx = lx, y, ly = ly, prog = FALSE)$frequency
  k_j <- k_j^q
  nck_j <- length(k_j)

  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)$frequency
  k <- k^q
  nck <- length(k)

  # Transfer entropy
  #------------------------------
  fraction <- matrix(NA, nck1_j, 2)
  numerator <- matrix(NA, nck1, 1)
  denominator <- matrix(NA, nck1_j, 1)

  for(i in 1:nck1) {
    numerator[i, 1] <- k1[i] / sum(k)
  }

  for(i in 1:nck1_j){
    denominator[i, 1] <- k1_j[i] / sum(k_j)
  }


  ren_entropy <- 1 / (1 - q) *
    log2(colSums(numerator)[1] / colSums(denominator)[1])

  return(list(transentropy = ren_entropy,
              numclassk1_j = nck1_j,
              numclassk1 = nck1,
              numclassk_j = nck_j,
              numclassk = nck))
}
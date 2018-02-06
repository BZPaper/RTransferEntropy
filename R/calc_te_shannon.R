#' Function that calculates the transfer entropy between two time series x and
#' y. The information flow from y to x is measured. Change x and y in function
#' call to infer the dominant direction of the information flow. Calculated
#' transfer entropy measure is Shannon transfer entropy.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#' @param lx Markov order of x
#' @param ly Markov order of y
#'
#' @return returns a list
#' @keywords internal
#' @export
#'
#' @examples
#'
calc_te_shannon <- function(x,
                            lx,
                            y,
                            ly) {

  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, lx = lx, y, ly = ly)$frequency
  nck1_j <- length(k1_j)

  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)$frequency
  nck1 <- length(k1)

  # x(k) and y(j)
  k_j <- cluster_gen(x, lx = lx, y, ly = ly, prog = FALSE)$frequency
  nck_j <- length(k_j)

  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)$frequency
  nck <- length(k)

  # Transfer entropy
  #------------------------------
  entropy <- numeric(nck1_j)
  for(i in 1:nck1_j){
    p1 <- k1[paste(strsplit(names(k1_j[i]), " ")[[1]][1:(lx + 1)],
                   collapse = " ")]
    p2 <- k_j[paste(strsplit(names(k1_j[i]), " ")[[1]][-(lx + 1)],
                    collapse = " ")]
    p3 <- k[paste(strsplit(names(k1_j[i]), " ")[[1]][1:lx],
                  collapse = " ")]
    entropy[i] <- k1_j[i] * log2((k1_j[i] * p3) / (p2 * p1))
  }

  return(list(transentropy = sum(entropy),
              numclassk1_j = nck1_j,
              numclassk1 = nck1,
              numclassk_j = nck_j,
              numclassk = nck))
}

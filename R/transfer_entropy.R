#' Function that calculates the transfer entropy between two time series x and 
#' y. The information flow from y to x is measured. Change x and y in function 
#' call to infer the dominant direction of the information flow.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#' @param lx x(t)
#' @param ly y(t)
#'
#' @return returns a list
#' @export
#'
#' @examples
#' 
transfer_entropy <- function(x, y, lx, ly){
  
  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, y, lx = lx, ly = ly)$frequency
  nck1_j <- length(k1_j)
  
  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)$frequency
  nck1 <- length(k1)
  
  # x(k) and y(j)
  k_j <- cluster_gen(x, y, lx = lx, ly = ly, prog = FALSE)$frequency
  nck_j <- length(k_j)
  
  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)$frequency
  nck <- length(k)
  
  # Transfer entropy
  entropy <- numeric(nck1_j)
  for(i in 1:nck1_j){
    p1 <- k1[substr(names(k1_j[i]), 1, (lx + 1))]
    p2 <- k_j[paste(unlist(strsplit(names(k1_j[i]), split = NULL))[-(lx + 1)],
                    collapse = "")]
    p3 <- k[substr(names(k1_j[i]), 1, lx)]
    entropy[i] <- k1_j[i] * log2((k1_j[i] * p3) / (p2 * p1))
  }
  
  return(list(transentropy = sum(entropy),
              numclassk1_j = nck1_j,
              numclassk1 = nck1,
              numclassk_j = nck_j,
              numclassk = nck))
}
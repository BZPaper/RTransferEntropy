#' Function for bootstrapping transfer entropy under H1.
#'
#' @param x a vector of coded values
#' @param valuestab table of codes and associated values
#' @param shuffle if TRUE, shuffled transfer entropy is calculated
#' @param lx x(k)
#' @param ly y(j)
#' @param const if TRUE, then shuffle is constant for all bootstraps
#' @param constx 
#' @param consty 
#' @param nrep 
#' @param shuffles 
#'
#' @return returns a vector
#' @export
#'
#' @examples
#' 
trans_boot_H1 <- function(x,
                          valuestab, 
                          shuffle = TRUE, 
                          lx, 
                          ly, 
                          const = FALSE, 
                          constx = NULL, 
                          consty = NULL, 
                          nrep = 10,
                          shuffles = 6) {
  
  # Bootstrap
  bootx <-  Markov_boot_step(x, lx) 
  
  # Decoding
  dmat <- decode_sample(bootx, valuestab)
  
  #Transfer entropy of bootstrapped process (with shuffling)
  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- transfer_entropy(dmat[, 1], dmat[, 2], lx = lx, 
                               ly = ly)$transentropy - constx
      # Lead = y
      dtey <- transfer_entropy(dmat[, 2], dmat[, 1], lx = ly, 
                               ly = lx)$transentropy - consty
    } else {
      # Lead = x
      dtex <- shuffled_transfer_entropy(nrep, 
                                        shuffles,
                                        diff = TRUE,
                                        dmat[, 1],
                                        lx = lx, 
                                        dmat[, 2],
                                        ly = ly)
      # Lead = y
      dtey <- shuffled_transfer_entropy(nrep, 
                                        shuffles,
                                        diff = TRUE,
                                        dmat[, 2],
                                        lx = ly, 
                                        dmat[, 1],
                                        ly = lx)
    }
  } else {
    # Transfer entropy of bootstrapped process (without shuffling)
    # Lead = x
    dtex <- transfer_entropy(dmat[, 1], dmat[, 2], lx = lx, 
                             ly = ly)$transentropy
    # Lead = y
    dtey <- transfer_entropy(dmat[, 2], dmat[, 1], lx = ly, 
                             ly = lx)$transentropy
  }
  
  teboot <- c(dtex, dtey) 
  names(teboot) <- c("dtex", "dtey")
  
  return(teboot)
}
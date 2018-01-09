#' Function for bootstrapping Rényi transfer entropy under H0.
#'
#' @param x a vector of coded values
#' @param lx x(k)
#' @param y a vector of coded values
#' @param ly y(j)
#' @param q weighting parameter in Rényi transfer entropy
#' @param shuffle 
#' @param const if TRUE, then shuffle is constant for all bootstraps
#' @param constx 
#' @param consty 
#' @param nreps 
#' @param shuffles 
#' @param ncores 
#'
#' @return returns a vector
#' @export
#'
#' @examples
#' 
trans_boot_H0_ren <- function(x,
                              lx,
                              y,
                              ly,
                              q,
                              shuffle = TRUE,
                              const = FALSE,
                              constx = NULL,
                              consty = NULL,
                              nreps = 2,
                              shuffles = 6,
                              ncores = parallel::detectCores() - 1) {
  
  bootx <- Markov_boot_step(x, lx)
  booty <- Markov_boot_step(y, ly)
  
  if (shuffle) {
    if (const) {
      # Lead = x
      dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy 
      - constx
      # Lead = y
      dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy 
      - consty
    } else {
      constx <- shuffled_transfer_entropy_ren(bootx,
                                              lx = lx,
                                              y,
                                              ly = ly,
                                              q,
                                              nreps,
                                              shuffles,
                                              diff = TRUE,
                                              ncores)
      consty <- shuffled_transfer_entropy_ren(booty,
                                              lx = ly,
                                              x,
                                              ly = lx,
                                              q,
                                              nreps,
                                              shuffles,
                                              diff = TRUE,
                                              ncores)
      
      # Lead = x
      dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy 
      - constx
      # Lead = y
      dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy 
      - consty
    }
  } else {
    # Lead = x
    dtex <- transfer_entropy_ren(bootx, lx = lx, y, ly = ly, q)$transentropy
    # Lead = y
    dtey <- transfer_entropy_ren(booty, lx = ly, x, ly = lx, q)$transentropy
  }
  
  teboot <- c(dtex, dtey)
  names(teboot) <- c("dtex", "dtey")
  
  return(teboot)
}
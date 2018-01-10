#' Function to implement Rényi transfer entropy.
#'
#' @param x 
#' @param lx1 
#' @param y 
#' @param ly 
#' @param q 
#' @param shuffle 
#' @param const 
#' @param nreps 
#' @param shuffles 
#' @param ncores 
#' @param quantiles 
#' @param bins 
#' @param limits 
#' @param boots 
#' @param nboot 
#'
#' @return returns a list
#' @export
#'
#' @examples
#' 
Renyi_transfer_entopy <- function(x,
                                  lx1,
                                  y,
                                  ly,
                                  q,
                                  shuffle = TRUE,
                                  const = FALSE,
                                  nreps = 2,
                                  shuffles = 6,
                                  ncores = parallel::detectCores() - 1,
                                  quantiles = c(5, 95),
                                  bins = NULL,
                                  limits = NULL,
                                  boots,
                                  nboot){
  
  # Code time series
  x <- code_sample(x, type = "quantiles", quantiles, bins, limits)
  y <- code_sample(y, type = "quantiles", quantiles, bins, limits)
  
  # Calculate transfer entropy (withour shuffling)
  # Lead = x
  tex <- transfer_entropy_ren(x, lx = lx, y, ly = ly, q)$transentropy
  # Lead = y
  tey <- transfer_entropy_ren(y, lx = ly, x, ly = lx, q)$transentropy
  
  # Calculate transfer entropy (with shuffling)
  constx <- shuffled_transfer_entropy_ren(x,
                                          lx = lx,
                                          y,
                                          ly = ly,
                                          q,
                                          nreps,
                                          shuffles,
                                          diff = FALSE,
                                          ncores)
  
  consty <- shuffled_transfer_entropy_ren(y,
                                          lx = ly,
                                          x,
                                          ly = lx,
                                          q,
                                          nreps,
                                          shuffles,
                                          diff = FALSE,
                                          ncores)
  
  # Lead = x
  stex <- tex - constx
  # Lead = y
  stey <- tey - consty
  
  # Bootstrap
  cl <- parallel::makeCluster(ncores)
  on.exit({
    parallel::stopCluster(cl)
  })
  
  parallel::clusterExport(cl, c("x", "y", "n", "lx", "ly", "q", "nreps", 
                                "nboot"),
                          envir = environment())
  
  seeds <- rnorm(boots)
  
  boot <- parallel::parLapply(cl, seeds, function(seed) {
    set.seed(seed)
    res <- replicate(nboot,
                     trans_boot_H0_ren(x,
                                       lx, 
                                       y, 
                                       ly, 
                                       q, 
                                       shuffle = TRUE,
                                       const = TRUE,
                                       constx = 0,
                                       consty = 0,
                                       nreps))
    return(res)
  })
  
  ste <- mean(unlist(boot))
  
  return(list(tex   = tex,
              tey   = tey,
              S_tex = S_tex,
              S_tey = S_tey,
              bootstrap_H0 = boot))
}
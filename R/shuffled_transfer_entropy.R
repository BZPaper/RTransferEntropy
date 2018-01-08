#' Function to calculate the effective transfer entropy, as given by the
#' difference between the transfer entropy calculated from a sample and the
#' respective shuffled transfer entropy.
#'

#' @param x a vector of coded values
#' @param lx x(k)
#' @param y a vector of coded values
#' @param ly y(j)
#' @param ncores the number of cores
#' @param nreps number of replications per shuffle
#' @param shuffles number of shuffles
#' @param diff if TRUE, the effective transfer entropy is calculated, otherwise
#' only the shuffled transfer entropy is returned
#'
#' @return returns a numeric scalar
#' @export
#'
#' @examples
#'
shuffled_transfer_entropy <- function(x,
                                      lx,
                                      y,
                                      ly,
                                      nreps = 2,
                                      shuffles = 6,
                                      diff = TRUE,
                                      ncores = parallel::detectCores() - 1) {

  n <- length(x)

  cl <- parallel::makeCluster(ncores)
  on.exit({
    parallel::stopCluster(cl)
  })

  parallel::clusterExport(cl, c("nreps", "x", "y", "n", "lx", "ly"),
                          envir = environment())

  seeds <- rnorm(shuffles)

  shuffle <- parallel::parLapply(cl, seeds, function(seed) {
    set.seed(seed)
    res <- replicate(nreps,
                     transfer_entropy(x = x,
                                      y = sample(y, n, replace = TRUE),
                                      lx = lx, ly = ly)$transentropy)
    return(res)
  })

  ste <- mean(unlist(shuffle))

  if (diff) {
    te <- transfer_entropy(x = x, y = y, lx = lx, ly = ly)$transentropy - ste
  } else {
    te <- ste
  }

  return(te)
}

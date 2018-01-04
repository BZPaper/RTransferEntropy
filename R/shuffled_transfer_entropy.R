#' Function to calculate the effective transfer entropy, as given by the
#' difference between the transfer entropy calculated from a sample and the
#' respective shuffled transfer entropy.
#'
#' @param nrep
#' @param shuffles
#' @param diff if TRUE, the effective transfer entropy is calculated, otherwise
#' only the shuffled transfer entropy is returned
#' @param x a vector of coded values
#' @param lx x(k)
#' @param y a vector of coded values
#' @param ly y(j)
#' @param ncores the number of cores
#'
#' @return returns a numeric scalar
#' @export
#'
#' @examples
#'
shuffled_transfer_entropy <- function(nrep = 10,
                                      shuffles = 6,
                                      diff = TRUE,
                                      x,
                                      lx,
                                      y,
                                      ly,
                                      ncores = parallel::detectCores() - 1) {

  n <- length(x)
  nreps <- round((nrep + 1) / shuffles)

  cl <- parallel::makeCluster(ncores)
  on.exit({
    parallel::stopCluster(cl)
  })

  shuffle <- parallel::parLapply(cl, 1:shuffles, function(i) {
    res <- replicate(nreps,
                     transfer_entropy(x = x,
                                      y = sample(y, n, replace = TRUE),
                                      lx = lx, ly = ly)$transentropy)
    return(res)
  })

  ste <- mean(as.numeric(as.character(unlist(shuffle))))

  if (diff) {
    te <- transfer_entropy(x = x, y = y, lx = lx, ly = ly)$transentropy - ste
  } else {
    te <- ste
  }

  return(te)
}

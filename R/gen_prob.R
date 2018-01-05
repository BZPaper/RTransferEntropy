#' Function generates probabilities for step-wise markov bootstrap.
#'
#' @param x a vector of coded values
#' @param lx x(k)
#'
#' @return returns a list
#' @export
#'
#' @examples
#'
gen_prob <- function(x,
                     lx) {

  n <- length(x)

  # Without lags
  px <- table(x)

  # With lags
  nclust <- n - lx
  clustlist <- list()

  for (i in 1:nclust) {
    clustlist[[i]] <- x[i:(i + lx)]
  }

  clust <- unlist(lapply(clustlist, function(x) paste(x, collapse = "")))
  mprob <- table(clust)
  nprob <- length(mprob)

  values <- names(mprob)
  pvec <- numeric(nprob)
  names(pvec) <- values

  for (i in 1:nprob) {
    p1 <- mprob[i]
    p2 <- px[substr(values[i], 1, 1)]
    pvec[i] <- p1/p2
  }

  transprob <- split(pvec, substr(values, 1, 1))

  return(list(px = px, transprob = transprob))
}

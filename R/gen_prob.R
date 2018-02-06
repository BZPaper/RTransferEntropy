#' Function generates state and transition probabilities for
#' step-wise markov bootstrap.
#'
#' @param x a vector of coded values
#' @param lx Markov order of x
#'
#' @return returns a list
#' @export
#' @keywords internal
#'
#' @examples
#'
gen_prob <- function(x,
                     lx) {

  n <- length(x)

  # Calculate state probabilities
  px <- table(x)/n

  # Calculate transition probabilities
  nclust <- n - lx
  clustlist <- list()

  for (i in 1:nclust) {
    clustlist[[i]] <- x[i:(i + lx)]
  }

  clust <- unlist(lapply(clustlist, function(x) paste(x, collapse = "")))
  mprob <- table(clust)/nclust
  nprob <- length(mprob)

  values <- names(mprob)
  pvec <- numeric(nprob)
  names(pvec) <- values

  for (i in 1:nprob) {
    p1 <- mprob[i]
    select <- substr(values[i], 1, 1)
    p2 <- sum(mprob[grep(paste("^", select, sep = ""), values, value = TRUE)])
    pvec[i] <- p1/p2
  }

  # Collect transition probabilities
  transprob <- split(pvec, substr(values, 1, 1))

  return(list(px = px, transprob = transprob))
}

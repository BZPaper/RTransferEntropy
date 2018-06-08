# Function calculates state and transition probabilities for step-wise
# bootstrap of the Markov chain. Bootstrapping the Markov chain is only
# relevant if statistical inference is required. State and transition
# probabilities are calculated based on relative frequencies and transition
# probabilities sum up to one. Function used internally.
#
gen_prob <- function(x,
                     lx) {

  n <- length(x)

  # Calculate state probabilities (relative frequencies of coded values)
  px <- table(x) / n

  # Calculate transition probabilities (relative frequencies of blocks
  # containing current and lagged values of coded time series)
  nclust <- n - lx
  clustlist <- list()

  for (i in 1:nclust) {
    clustlist[[i]] <- x[i:(i + lx)]
  }

  clust <- unlist(lapply(clustlist, function(x) paste(x, collapse = "")))
  mprob <- table(clust) / nclust
  nprob <- length(mprob)

  values <- names(mprob)
  pvec <- numeric(nprob)
  names(pvec) <- values

  for (i in 1:nprob) {
    p1 <- mprob[i]
    select <- substr(values[i], 1, 1)
    p2 <- sum(mprob[grep(paste("^", select, sep = ""), values, value = TRUE)])
    pvec[i] <- p1 / p2
  }

  # Collect transition probabilities
  transprob <- split(pvec, substr(values, 1, 1))

  return(list(px = px, transprob = transprob))
}

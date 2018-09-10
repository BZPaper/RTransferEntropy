# Used internally for the bootstrap of the Markov chain.
#
combine_sample <- function(x, y) {
  n <- length(x)
  newsample <- cbind(x, y)
  newsample <- apply(newsample, 1, function(x) paste(x, collapse = ""))

  values <- sort(unique(newsample))
  numvalues <- length(values)

  pool <- LETTERS
  valuestab <- pool[1:numvalues]
  names(valuestab) <- values

  csample <- numeric(n)
  csample[] <- NA

  if (length(pool) < numvalues) {
    print("Too many combinations.")
    stop
  } else {
    for (i in 1:numvalues) {
      csample <- ifelse(newsample == values[i], pool[i], csample)
    }
  }

  return(list(valuestab = valuestab, csample = csample))
}

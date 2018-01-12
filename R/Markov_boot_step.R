#' Function for step-wise bootstraping the Markov chain.
#'
#' @param x a vector of coded values
#' @param lx x(k)
#' @param burn number of observations that are dropped from the beginning of
#' the bootstrapped time series
#'
#' @return returns a vector of strings
#' @export
#'
#' @examples
#'
Markov_boot_step <- function (x,
                              lx,
                              burn = 50) {

  n <- length(x) + burn
  bootlist <- list()

  # First draw
  probs <- gen_prob(x, lx)
  set.seed(n)
  fdraw <- sample(names(probs$px), 1, prob = probs$px)
  bootlist[[1]] <- fdraw

  # Draws

  for (i in 2:(lx + 1)) {
    probs <- gen_prob(x, lx = i - 1)
    tprob <- probs$transprob
    set.seed(n)
    draw <- ifelse(length(names(tprob[[bootlist[[i - 1]]]])) == 1,
                   names(tprob[[bootlist[[i - 1]]]]),
                   sample(names(tprob[[bootlist[[i - 1]]]]), 1,
                          prob = tprob[[bootlist[[i - 1]]]]))
    bootlist[[i]] <- substr(draw, i, i)
  }

  tprob <- gen_prob(x, lx)$transprob

  for (i in (lx + 2):n) {
    set.seed(n)
    draw <- ifelse(length(names(tprob[[bootlist[[i - 1]]]])) == 1,
                   names(tprob[[bootlist[[i - 1]]]]),
                   sample(names(tprob[[bootlist[[i - 1]]]]), 1,
                          prob = tprob[[bootlist[[i - 1]]]]))
    bootlist[[i]] <- substr(draw, lx + 1, lx + 1)
  }

  bootvec <- unlist(bootlist)
  bootvec <- bootvec[(burn + 1):n]

  return(bootvec)
}

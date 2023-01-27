# Function used for sampling from the coded time series. New time series is
# used for the H0 bootstrap. Used internally by transfer_entropy; same arguments.
#
markov_boot_step <- function(x, lx, burn = 50) {
  n <- length(x) + burn
  bootvec <- numeric(n)

  # First draw
  pr <- freq_table(x[-(length(x) - lx + 1):-length(x)])
  lb <- sample(names(pr), 1, prob = pr)
  bootvec[1] <- lb

  # Draws
  for (i in 2:(lx + 1)) {
    tprob <- calculate_transition_probabilities(x, lx = i - 1)
    val <- tprob[[lb]]
    draw <- sample(names(val), 1, prob = val)
    lb <- strsplit(draw, " ")[[1]][i]
    bootvec[i] <- lb
  }

  tprob <- calculate_transition_probabilities(x, lx)

  for (i in (lx + 2):n) {
    val2 <- tprob[[lb]]

    # guard against an edge case
    if (is.null(val2) || val2[[1]] == 0) {
      # unlikely, therefore not likely / important
      draw <- sample(names(val2), 1)
      lb <- setdiff(strsplit(draw, " ")[[1]], lb)[1]
      val2 <- tprob[[lb]]
    }
    val <- val2

    if (length(val2) < 1) print(str(val2))
    draw <- names(val2)[sample.int(length(val2), 1, prob = val2)]
    lb <- strsplit(draw, " ")[[1]][lx + 1]
    bootvec[i] <- lb
  }

  bootvec <- bootvec[(burn + 1):n]

  return(as.numeric(bootvec))
}

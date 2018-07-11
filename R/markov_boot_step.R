# Function used for sampling from the coded time series. New time series is
# used for the H0 bootstrap. Used internally by transfer_entropy; same arguments.
#
markov_boot_step <- function(x, lx, burn = 50) {

  n <- length(x) + burn
  bootvec <- numeric(n)

  # First draw
  pr <- freq_table(x)
  bootvec[1] <- sample(names(pr), 1, prob = pr)

  # Draws
  for (i in 2:(lx + 1)) {
    tprob <- calculate_transition_probabilities(x, lx = i - 1)
    val <- tprob[[bootvec[i - 1]]]

    draw <- ifelse(length(names(val)) == 1,
                   names(val),
                   sample(names(val), 1, prob = val))

    bootvec[i] <- strsplit(draw, " ")[[1]][i]
  }

  tprob <- calculate_transition_probabilities(x, lx)

  for (i in (lx + 2):n) {
    val <- tprob[[bootvec[i - 1]]]

    draw <- ifelse(length(names(val)) == 1,
                   names(val),
                   sample(names(val), 1, prob = val))

    bootvec[i] <- strsplit(draw, " ")[[1]][lx + 1]
  }

  bootvec <- bootvec[(burn + 1):n]

  return(as.numeric(bootvec))
}

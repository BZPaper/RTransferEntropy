# Function calculates state and transition probabilities for step-wise
# bootstrap of the Markov chain. Bootstrapping the Markov chain is only
# relevant if statistical inference is required. State and transition
# probabilities are calculated based on relative frequencies and transition
# probabilities sum up to one. Function used internally.
#
gen_prob <- function(x, lx) {

  if (is.numeric(x)) x <- as.character(x)

  px <- freq_table(x)
  transprob <- calculate_transition_probabilities(x, lx)

  return(list(px = px, transprob = transprob))
}

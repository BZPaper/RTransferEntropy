#' Function to implement Shannon transfer entropy.
#'
#' @param x a vector of coded values
#' @param lx Markov order of x
#' @param y a vector of coded values
#' @param ly Markov order of y
#' @param const if TRUE, then shuffle is constant for all bootstraps
#' @param constx constant value substracted from transfer entropy measure
#' @param consty constant value substracted from transfer entropy measure
#' @param shuffle if TRUE, shuffled transfer entropy is calculated
#' @param nreps number of replications for each shuffle
#' @param shuffles number of shuffles
#' @param ncores number of cores in parallel computation
#' @param type bins, limits or quantiles of empirical distribution to discretize
#' the data
#' @param quantiles quantiles to use for discretization
#' @param bins the number of bins with equal width used for discretization
#' @param limits limits used for discretization
#' @param parcalc "yes" if parallel computations for the Markov chain bootstrap
#' @param boots number of bootstrap samples
#' @param nboot number of bootstrap replications
#' @param burn number of observations that are dropped from the beginning of
#' the bootstrapped Markov chain
#'
#' @return returns a list
#' @export
#'
#' @examples
#'
Shannon_transfer_entopy <- function(x,
                                    lx,
                                    y,
                                    ly,
                                    shuffle = TRUE,
                                    const = FALSE,
                                    constx = 0,
                                    consty = 0,
                                    nreps = 2,
                                    shuffles = 6,
                                    ncores = parallel::detectCores() - 1,
                                    type = "quantiles",
                                    quantiles = c(5, 95),
                                    bins = NULL,
                                    limits = NULL,
                                    parcalc = "yes",
                                    boots = 1,
                                    nboot,
                                    burn = 50) {

  # Code time series
  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # Calculate transfer entropy (without shuffling)
  # Lead = x
  tex <- transfer_entropy(x, lx = lx, y, ly = ly)$transentropy
  # Lead = y
  tey <- transfer_entropy(y, lx = ly, x, ly = lx)$transentropy

  # Calculate transfer entropy (with shuffling)
  constx <- shuffled_transfer_entropy(x,
                                      lx = lx,
                                      y,
                                      ly = ly,
                                      nreps,
                                      shuffles,
                                      diff = FALSE,
                                      ncores)

  consty <- shuffled_transfer_entropy(y,
                                      lx = ly,
                                      x,
                                      ly = lx,
                                      nreps,
                                      shuffles,
                                      diff = FALSE,
                                      ncores)

  # Lead = x
  stex <- tex - constx
  # Lead = y
  stey <- tey - consty

  # Bootstrap
  if (parcalc == "yes") {
    cl <- parallel::makeCluster(ncores)
    on.exit({
      parallel::stopCluster(cl)
    })

    parallel::clusterExport(cl, c("nreps", "x", "y", "n", "lx", "ly",
                                  "code_sample", "cluster_gen",
                                  "transfer_entropy",
                                  "shuffled_transfer_entropy", "gen_prob",
                                  "Markov_boot_step", "trans_boot_H0"),
                            envir = environment())

    seeds <- rnorm(boots)

    boot <- parallel::parLapply(cl, seeds, function(seed) {
      set.seed(seed)
      res <- replicate(nboot,
                       trans_boot_H0(x,
                                     lx = lx,
                                     y,
                                     ly = ly,
                                     burn,
                                     shuffle,
                                     const,
                                     constx,
                                     consty,
                                     nreps,
                                     shuffles,
                                     ncores))
      return(res)
    })
  } else {
    boot <- replicate(nboot,
                      trans_boot_H0(x,
                                    lx = lx,
                                    y,
                                    ly = ly,
                                    burn,
                                    shuffle,
                                    const,
                                    constx,
                                    consty,
                                    nreps,
                                    shuffles,
                                    ncores))
  }


  return(list(tex   = tex,
              tey   = tey,
              stex = stex,
              stey = stey,
              bootstrap_H0 = boot))
}

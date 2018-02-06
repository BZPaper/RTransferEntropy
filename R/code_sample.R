#' Codes the sample by assigning data to bins that are, for example, based on
#' the quantiles of the empirical distribution of the sample.
#'
#' @param x a vector of numerical values
#' @param type bins, limits or quantiles of empirical distribution to discretize
#' the data
#' @param quantiles quantiles to use for discretization
#' @param bins the number of bins with equal width used for discretization
#' @param limits limits used for discretization
#' @param scale a scale parameter
#'
#' @return returns a numerical vector
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#'  set.seed(42)
#'  x <- rnorm(100)
#'  code_sample(x)
#'
#'  # or coming from a data.table framework:
#'  library(data.table)
#'  set.seed(42)
#'  dt <- data.table(x = rnorm(100))
#'  dt[, sample := code_sample(dt$x)]
#' }
code_sample <- function(x,
                        type = "quantiles",
                        quantiles = c(5, 95),
                        bins = NULL,
                        limits = NULL,
                        scale = 1e10) {

  if (type %in% c("bins", "limits")) {
    UB <- max(x)
    LB <- min(x)

    # find the respective OSeq for the time series
    if (type == "bins") {
      if(is.null(bins)) {
        stop(cat(paste("Warning: Bins not defined", "\n",
                       "Execution halted","\n")))
      }
      OSeq <- LB + ((UB - LB) / bins) * (0:(bins))
    } else {
      if(is.null(limits)) {
        stop(cat(paste("Warning: Limits not defined", "\n",
                       "Execution halted","\n")))
      }
      limits <- sort(limits)
      OSeq <- c(LB, limits, UB)
    }

    OSeq[length(OSeq)] <- UB + 1

    for (j in 1:(length(OSeq) - 1)) {
      x[x >= OSeq[j] & x < OSeq[j + 1]] <- j * scale
    }
  } else if (type == "quantiles") {
    Qtl <- quantile(x, type = 8, probs = quantiles / 100)
    Qtl <- c(Qtl, max(x))
    Qlength <- length(Qtl)

    for (j in 1:Qlength){
      x[x <= Qtl[j]] <- j * scale
    }
  }
  x <- x / scale

  return(x)
}

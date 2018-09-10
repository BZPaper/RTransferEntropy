# Codes the time series by discretizing the data, i.e. assigning values to bins
# that are either based on the quantiles of the empirical distribution of the
# given time series or are specified directly over a chosen number of bins or
# predefined value-limits. Function used internally.
#
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
      if (is.null(bins)) stop("Warning: Bins not defined\nExecution halted\n")

      OSeq <- LB + ((UB - LB) / bins) * (0:(bins))
    } else {
      if (is.null(limits)) stop("Warning: Bins not defined\nExecution halted\n")

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

    for (j in 1:Qlength) {
      x[x <= Qtl[j]] <- j * scale
    }
  }
  x <- x / scale

  return(x)
}

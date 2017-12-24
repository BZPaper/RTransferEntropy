#' Codes a vector into different bins
#'
#' Codes the sample by assigning data to bins that are, for example, based on
#' the quantiles of the empirical distribution of the sample.
#'
#' @param Data Column of a data table
#' @param Type Assign a certain number of bins, limits or use quantiles of
#'             empirical distribution to tansform data into discrete form
#' @param Quantile Define the quantiles to use for discretization
#' @param Bins Define the number of bins with equal width used for
#'             discretization
#' @param Limits Define limits used for discretization explicitly
#' @param Scale a scale parameter
#'
#' @return Function returns a data table column
#' @export
#'
#' @examples
#'
code_sample <- function(Data,
                        Type = "quantiles",
                        Quantile = c(5, 95),
                        Bins = NULL,
                        Limits = NULL,
                        Scale = 1e10) {

  CodeData <- copy(Data)
  nam <- copy(colnames(Data))
  setnames(CodeData, nam, "TS")

  if (Type %in% c("bins", "limits")) {
    UB <- max(CodeData)
    LB <- min(CodeData)

    # find the respective OSeq for the time series
    if (Type == "bins") {
      OSeq <- LB + ((UB - LB) / Bins) * (0:(Bins))
    } else {
      Limits <- sort(Limits)
      OSeq <- c(LB, Limits, UB)
    }

    OSeq[length(OSeq)] <- UB + 1

    for (j in 1:(length(OSeq) - 1)) {
      CodeData <- CodeData[TS >= OSeq[j] & TS < OSeq[j + 1], TS := j * Scale]
    }
  } else if (Type == "quantiles") {
    Qtl <- quantile(CodeData[, TS], type = 8, probs = Quantile/100)
    Qtl <- c(Qtl, max(CodeData))
    Qlength <- length(Qtl)

    for (j in 1:Qlength){
      CodeData <- CodeData[TS <= Qtl[j], TS := j * Scale]
    }
  }

  CodeData <- CodeData[, TS := TS / Scale]
  CodeData <- CodeData[order(TS)]
  setnames(CodeData, "TS", nam)

  return(CodeData)
}

#' @useDynLib RTransferEntropy
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile sd
NULL

star <- function(x) {
  ifelse(is.null(x) || is.na(x), "",
         ifelse(x < 0.001, "***",
                ifelse(x < 0.01, "**",
                       ifelse(x < 0.05, "*",
                              ifelse(x < 0.1, ".", "")))))
}

fupper <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

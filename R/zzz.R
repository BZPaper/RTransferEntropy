#' @useDynLib RTransferEntropy
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile sd
#' @importFrom future plan multisession sequential
NULL

.onAttach <- function(...) {
  set_quiet(FALSE)
}

# for some p-values (x) return the stars
star <- function(x) {
  ifelse(is.null(x) || is.na(x), "",
         ifelse(x < 0.001, "***",
                ifelse(x < 0.01, "**",
                       ifelse(x < 0.05, "*",
                              ifelse(x < 0.1, ".", "")))))
}

# first to upper
fupper <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

# returns an s if n > 1 (i.e., sprintf("we have n = %s sample%s", n, mult_s(n)))
mult_s <- function(n) ifelse(n > 1, "s", "")

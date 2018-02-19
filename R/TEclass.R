#' Prints a transfer-entropy result
#'
#' @param x a TEResult
#' @param ...
#'
#' @return invisible the text
#' @export
#'
#' @examples
#' # see ?transfer_entropy()
print.TEResult <- function(x, ...) {
  star <- function(x) {
    ifelse(is.null(x) || is.na(x), "",
           ifelse(x < 0.001, "***",
                  ifelse(x < 0.01, "**",
                         ifelse(x < 0.05, "*", ""))))
  }
  n <- 64
  str <- c(
    "Transfer Entropy Result:",
    sprintf("Direction  %10s  %10s  %10s  %10s  %5s",
            "te", "ete", "se", "p-value", "sig"),
    paste(rep("-", n), collapse = ""),
    sprintf("X->Y       %10.05f  %10.05f  %10.05f  %10.05f  %5s",
            x$te_xy, x$ete_xy, x$se_xy, x$p_xy, star(x$p_xy)),
    sprintf("Y->X       %10.05f  %10.05f  %10.05f  %10.05f  %5s",
            x$te_yx, x$ete_yx, x$se_yx, x$p_yx, star(x$p_yx)),
    paste(rep("-", n), collapse = ""),
    sprintf("Number of Observations: %s", x$nobs),
    paste(rep("-", n), collapse = "")
  )
  str <- paste(str, collapse = "\n")
  cat(str)
  return(invisible(str))
}

#' Prints a summary of a transfer-entropy result
#'
#' @param x a TEResult
#' @param ... additional arguments
#'
#' @return invisible the text
#' @export
#'
#' @examples
#' # see ?transfer_entropy
summary.TEResult <- function(x, ...) {
  str <- c(
    "Transfer Entropy Result:",
    sprintf("Number of Observations: %s", length(x$nobs))
  )
  str <- paste(str, collapse = "\n")
  cat(str)
  return(invisible(str))
}

#' Checks if an object is a TEResult
#'
#' @param x an object
#'
#' @return a boolean value if x is a TEResult
#' @export
#'
#' @examples
#' # see ?transfer_entropy
is.TEResult <- function(x) {
  inherits(x, "TEResult")
}

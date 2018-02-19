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
print.TEResult <- function(x, digits = 4, ...) {
  star <- function(x) {
    ifelse(is.null(x) || is.na(x), "",
           ifelse(x < 0.001, "***",
                  ifelse(x < 0.01, "**",
                         ifelse(x < 0.05, "*",
                                ifelse(x < 0.1, ".", "")))))
  }
  fupper <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

  # the number of chars per reported value
  n_digits <- max(10, digits + 2)
  # the width (in chars) of the overall output
  n <- 9 + (2 + n_digits) * 4 + 3 +2
  line <- paste(rep("-", n), collapse = "")

  nr_format <- sprintf("  %%%s.0%sf", n_digits, digits)
  header_els <- paste0("%9s", paste(rep(sprintf("  %%%ss", n_digits), 4), collapse = ""), "  %s")
  row_els <- paste0("%9s", paste(rep(nr_format, 4), collapse = ""), "  %s")

  str <- c(
    paste(fupper(x$entropy), "Transfer Entropy Results:"),
    line,
    sprintf(header_els,
            "Direction", "te", "ete", "se", "p-value", "sig"),
    line,
    sprintf(row_els, "X->Y", x$te_xy, x$ete_xy, x$se_xy, x$p_xy, star(x$p_xy)),
    sprintf(row_els, "Y->X", x$te_yx, x$ete_yx, x$se_yx, x$p_yx, star(x$p_yx)),
    line,
    paste0(sprintf("Number of Observations: %s", x$nobs),
           ifelse(x$entropy == "renyi", sprintf("\nQ: %s", x$q), "")),
    line,
    "Sig. P-values:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1"
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

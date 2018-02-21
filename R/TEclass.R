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

  # the number of chars per reported value
  n_digits <- max(10, digits + 2)

  # the width (in chars) of the overall output
  # ncol(x) - 1 of n_digits + 2 spaces (stars excluded)
  # plus the names (n_digits + 2 spaces)
  # plus 5 for the stars
  n <- (n_digits + 2) * (ncol(x$coef) + 1) + 5

  line <- paste(rep("-", n), collapse = "")

  # crate the header
  header_lengths <- c(rep(10, ncol(x$coef) + 1), 5)
  header_names <- c("Direction", "TE", "Eff. TE", "Std.Err.", "p-value", "sig")
  header <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
                         l = header_lengths, t = header_names),
                  collapse = "  ")

  str <- c(
    paste(fupper(x$entropy), "Transfer Entropy Results:"),
    line,
    header,
    line,
    textify_coef(x$coef, digits, 10),
    line,
    paste0(sprintf("Number of Observations: %s", x$nobs),
           ifelse(x$entropy == "renyi", sprintf("\nQ: %s", x$q), "")),
    line,
    "p-values: < 0.001 ‘***’, < 0.01 ‘**’, < 0.05 ‘*’, < 0.1 ‘.’"
  )
  str <- paste(str, collapse = "\n")
  cat(str)
  return(invisible(str))
}

# mat the matrix that contains the coefficients
# n the number of digits for the coefficients
# w the width of each number-field, defaults to 10
textify_coef <- function(mat, n, w = 10) {

  w <- max(10, n + 2)
  nr_fmt <- sprintf("%%%s.%sf", w, n)
  txt_fmt <- sprintf("%%%ss", w)

  # for each row, for each col, paste the number in the right format and
  # add the stars at the end
  txt <- apply(mat, 1, function(row_el) {
    res <- sapply(row_el, function(x) sprintf(nr_fmt, x))

    paste(c(res, sprintf("%5s", star(row_el[length(row_el)]))),
          collapse = "  ")
  })

  paste(sprintf(txt_fmt, names(txt)), txt, sep = "  ")
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

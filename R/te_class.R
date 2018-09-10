#' Prints a transfer-entropy result
#'
#' @param x a transfer_entropy
#' @param digits the number of digits to display, defaults to 4
#' @param boot if the bootstrapped results should be printed, defaults to TRUE
#' @param ... additional arguments, currently not in use
#'
#' @return invisible the text
#' @export
#'
#' @examples
#' # see ?transfer_entropy()
print.transfer_entropy <- function(x, digits = 4, boot = TRUE, ...) {

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

  # create the bootstrapped output:
  if (!is.matrix(x$boot) || !boot) {
    boot_res <- c(
      line,
      "For calculation of standard errors and p-values set nboot > 0"
    )
  } else {
    quants <- t(apply(x$boot, 1, function(b) quantile(b)))
    rownames(quants) <- c("X->Y", "Y->X")

    quant_hdr_l <- c(rep(8, ncol(quants) + 1))
    quant_hdr_n <- c("Direction", "0%", "25%", "50%", "75%", "100%")
    quant_hdr <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
                              l = quant_hdr_l, t = quant_hdr_n),
                       collapse = "  ")

    boot_res <- c(
      line,
      sprintf("Bootstrapped TE Quantiles (%s replications):", ncol(x$boot)),
      line,
      quant_hdr,
      line,
      textify_mat(quants, digits = digits, width = 8, stars = FALSE)
    )
  }

  text <- c(
    paste(fupper(x$entropy), "Transfer Entropy Results:"),
    line,
    header,
    line,
    textify_mat(x$coef, digits, 10),
    boot_res,
    line,
    paste0(sprintf("Number of Observations: %s", x$nobs),
           ifelse(x$entropy == "renyi", sprintf("\nQ: %s", x$q), "")),
    line,
    "p-values: < 0.001 '***', < 0.01 '**', < 0.05 '*', < 0.1 '.'"
  )
  text <- paste(text, collapse = "\n")
  cat(text)
  return(invisible(text))
}

# mat the matrix that contains the coefficients
# n the number of digits for the coefficients
# w the width of each number-field, defaults to 10
# stars if the last row represents the p-values and we want to calc the ***
textify_mat <- function(mat, digits, width = 10, stars = TRUE) {

  width <- max(width, digits + 2)
  nr_fmt <- sprintf("%%%s.%sf", width, digits)
  txt_fmt <- sprintf("%%%ss", width)

  # for each row, for each col, paste the number in the right format and
  # add the stars at the end
  txt <- apply(mat, 1, function(row_el) {
    res <- sapply(row_el, function(x) sprintf(nr_fmt, x))

    if (stars) {
      paste(c(res, sprintf("%5s", star(row_el[length(row_el)]))),
            collapse = "  ")
    } else {
      paste(res, collapse = "  ")
    }

  })

  paste0(" ", paste(sprintf(txt_fmt, names(txt)), txt, sep = "  "))
}

#' Prints a summary of a transfer-entropy result
#'
#' @param object a transfer_entropy
#' @param digits the number of digits to display, defaults to 4
#' @param ... additional arguments, passed to \code{\link[stats]{printCoefmat}}
#'
#' @return invisible the object
#' @export
#'
#' @examples
#' # see ?transfer_entropy
summary.transfer_entropy <- function(object, digits = 4, ...) {
  cat(sprintf("%s's Transfer Entropy\n\n", fupper(object$entropy)))
  cat("Coefficients:\n")
  printCoefmat(coef(object), ...)

  if (!is.matrix(object$boot)) {
    boot_res <- c(NULL)
  } else {
    quants <- t(apply(object$boot, 1, function(b) quantile(b)))
    rownames(quants) <- c("X->Y", "Y->X")

    quant_hdr_l <- c(rep(8, ncol(quants) + 1))
    quant_hdr_n <- c("Direction", "0%", "25%", "50%", "75%", "100%")
    quant_hdr <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
                              l = quant_hdr_l, t = quant_hdr_n),
                       collapse = "  ")

    boot_res <- c(
      sprintf("\nBootstrapped TE Quantiles (%s replications):",
              ncol(object$boot)),
      quant_hdr,
      textify_mat(quants, digits = digits, width = 8, stars = FALSE)
    )
    boot_res <- paste(boot_res, collapse = "\n")
    cat(boot_res,"\n")
  }

  cat(sprintf("\nNumber of Observations: %i", object$nobs))

  return(invisible(object))
}

#' Checks if an object is a transfer_entropy
#'
#' @param x an object
#'
#' @return a boolean value if x is a transfer_entropy
#' @export
#'
#' @examples
#' # see ?transfer_entropy
is.transfer_entropy <- function(x) {
  inherits(x, "transfer_entropy")
}

#' Extract the Coefficient Matrix from a transfer_entropy
#'
#' @param object a transfer_entropy
#' @param ... additional arguments, currently not in use
#'
#' @return a Matrix containing the coefficients
#' @export
#'
#' @examples
#' set.seed(1234567890)
#' n <- 500
#' x <- rep(0, n + 1)
#' y <- rep(0, n + 1)
#'
#' for (i in seq(n)) {
#'   x[i + 1] <- 0.2 * x[i] + rnorm(1, 0, 2)
#'   y[i + 1] <- x[i] + rnorm(1, 0, 2)
#' }
#'
#' x <- x[-1]
#' y <- y[-1]
#'
#' te_result <- transfer_entropy(x, y, nboot = 100)
#' coef(te_result)
coef.transfer_entropy <- function(object, ...) {
  if (!is.transfer_entropy(object)) stop("object must be a transfer_entropy")
  return(object$coef)
}

# for some p-values (x) return the stars
star <- function(x) {
  ifelse(is.null(x) || is.na(x), "",
         ifelse(x < 0.001, "***",
                ifelse(x < 0.01, "**",
                       ifelse(x < 0.05, "*",
                              ifelse(x < 0.1, ".", "")))))
}

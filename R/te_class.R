#' Prints a transfer-entropy result
#'
#' @param x a transfer_entropy
#' @param digits the number of digits to display, defaults to 4
#' @param boot if the bootstrapped results should be printed, defaults to TRUE
#' @param probs numeric vector of quantiles for the bootstraps
#' @param ... additional arguments, currently not in use
#' @param tex if the data should be outputted as a TeX-string
#' @param ref the reference string of the LaTeX table (label) applies only if table = TRUE and tex = TRUE, defaults to FALSE
#' @param file a file where the results are printed to
#' @param table if the table environment should be printed as well (only applies if tex = TRUE), defaults to TRUE
#'
#' @return invisible the text
#' @export
#'
#' @examples
#' # construct two time-series
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
#' # Calculate Shannon's Transfer Entropy
#' te_result <- transfer_entropy(x, y, nboot = 100)
#'
#' print(te_result)
#'
#' # change the number of digits
#' print(te_result, digits = 10)
#'
#' # disable boot-print
#' print(te_result, boot = FALSE)
#'
#' # specify the quantiles of the bootstraps
#' print(te_result, probs = c(0, 0.1, 0.4, 0.5, 0.6, 0.9, 1))
#'
#' # get LaTeX output:
#' print(te_result, tex = TRUE)
#'
#' # set the reference label for LaTeX table
#' print(te_result, tex = TRUE, ref = "tab:te_result")
#'
#' \dontrun{
#' # file output
#' print(te_result, file = "te_result_file.txt")
#' print(te_result, tex = TRUE, file = "te_result_file.tex")
#' }
print.transfer_entropy <- function(x, digits = 4, boot = TRUE,
                                   probs = c(0, 0.25, 0.5, 0.75, 1),
                                   tex = FALSE,
                                   ref = NA,
                                   file = NA,
                                   table = TRUE,
                                   ...) {

  # the number of chars per reported value

  # crate the header
  val_names <- c("TE", "Eff. TE", "Std.Err.", "p-value", "sig")
  header_names <- c("Direction", val_names)

  n_digits <- max(max(nchar(val_names)), digits + 2)
  header_lengths <- c(9, rep(n_digits, ncol(x$coef)), 5)

  header <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
    l = header_lengths, t = header_names
  ),
  collapse = ifelse(tex, " & ", "  ")
  )

  line <- paste0(rep("-", max(nchar(header)), 59), collapse = "")
  # 59 chars in the p-value footnote

  # create the bootstrapped output:
  if (!is.matrix(x$boot)) {
    boot_res <- c(
      line,
      "For calculation of standard errors and p-values set nboot > 0"
    )
  } else if (!boot) {
    boot_res <- NULL
  } else {
    quants <- t(apply(x$boot, 1, function(b) quantile(b, probs = probs)))
    rownames(quants) <- c("X->Y", "Y->X")

    probs_nam <- paste0(probs * 100, "%")
    boot_hd_nam <- c("Direction", probs_nam)
    boot_hd_len <- c(
      9,
      rep(
        max(nchar(probs_nam), digits + 2),
        length(probs)
      )
    )

    boot_hd <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
      l = boot_hd_len, t = boot_hd_nam
    ),
    collapse = ifelse(tex, " & ", "  ")
    )

    line_width <- max(nchar(header), nchar(boot_hd), 59)
    # 59 chars in the p-value footnote

    line <- paste(rep("-", line_width), collapse = "")

    if (tex) {
      boot_res <- c(
        "  \\begin{subtable}[t]{\\linewidth}",
        "    \\centering",
        "    \\vspace{0pt}",
        "    \\caption{Panel B: Bootstrapped TE Quantiles}",
        "    \\begin{tabular}{rrrrrr}",
        "      \\toprule",
        paste("      ", gsub("%", "\\\\%", boot_hd), " \\\\"),
        "      \\midrule",
        paste("      ",
              textify_mat(quants, digits = digits, width = boot_hd_len, stars = FALSE, tex = tex)),
        "      \\bottomrule",
        sprintf("      \\multicolumn{6}{l}{Number of Replications: %d}", ncol(x$boot)),
        "    \\end{tabular}",
        "  \\end{subtable}"
      )

    } else {
      boot_res <- c(
        line,
        sprintf("Bootstrapped TE Quantiles (%s replications):", ncol(x$boot)),
        line,
        boot_hd,
        line,
        textify_mat(quants, digits = digits, width = boot_hd_len, stars = FALSE)
      )
    }
  }


  if (tex) {
    if (!is.na(ref) && ref == FALSE) ref <- NA

    text <- c(
      "% Make sure to include the following packages:",
      "% \\usepackage{booktabs}",
      "% \\usepackage{subcaption}",
      "",
      "% Table created by the RTransferEntropy Package",
      sprintf("%% Time of creation: %s", Sys.time()),
      "",
      ifelse(table, "\\begin{table}", ""),
      ifelse(table, sprintf("  \\caption{%s Transfer Entropy Result}", fupper(x$entropy)), ""),
      ifelse(!is.na(ref) && table, sprintf("  \\label{%s}", ref), ""),
      "  \\begin{subtable}[t]{\\linewidth}",
      "    \\centering",
      "    \\vspace{0pt}",
      "    \\caption{Panel A: Transfer Entropy Estimates}",
      "    \\begin{tabular}{rrrrrl}",
      "      \\toprule",
      paste("      ", header, " \\\\"),
      "      \\midrule",
      paste("      ", textify_mat(x$coef, digits, header_lengths, tex = tex)),
      "      \\bottomrule",
      ifelse(x$entropy == "shannon",
             "",
             sprintf("      \\multicolumn{6}{l}{Q: %s} \\\\", x$q)),
      sprintf("      \\multicolumn{6}{l}{Number of Observations: %d} \\\\", x$nobs),
      "      \\multicolumn{6}{l}{p-values: \\textless 0.001 '***', \\textless 0.01 '**', \\textless 0.05 '*', \\textless 0.1 '.' } \\\\",
      "    \\end{tabular}",
      "  \\end{subtable}",
      "",
      boot_res,
      ifelse(table, "\\end{table}", "")
    )

  } else {
    text <- c(
      paste(fupper(x$entropy), "Transfer Entropy Results:"),
      line,
      header,
      line,
      textify_mat(x$coef, digits, header_lengths),
      boot_res,
      line,
      paste0(
        sprintf("Number of Observations: %s", x$nobs),
        ifelse(x$entropy == "renyi", sprintf("\nQ: %s", x$q), "")
      ),
      line,
      "p-values: < 0.001 '***', < 0.01 '**', < 0.05 '*', < 0.1 '.'"
    )
  }

  text <- paste(text, collapse = "\n")

  if (is.na(file)) {
    cat(text, "\n")
  } else {
    writeLines(text, file)
    message("TE Result has been written to file: ", file)
  }

  return(invisible(text))
}

# mat the matrix that contains the coefficients
# n the number of digits for the coefficients
# w the width of each number-field, defaults to 10
# stars if the last row represents the p-values and we want to calc the ***
textify_mat <- function(mat, digits, width = 10, stars = TRUE, tex = FALSE) {

  # the first element is the direction (text)
  nr_fmt <- sprintf("%%%s.%sf", width[-1], digits)
  txt_fmt <- sprintf("%%%ss", width[1])

  if (stars) {
    if (ncol(mat) + 1 != length(nr_fmt)) {
      stop("width has to have the same lenghts as the num of columns of mat + 1")
    }
    star_fmt <- sprintf("%%%ss", width[length(width)])
    nr_fmt <- nr_fmt[-length(nr_fmt)]
  }

  clp <- ifelse(tex, " & ", "  ")
  # for each row, for each col, paste the number in the right format and
  # add the stars at the end
  txt <- apply(mat, 1, function(row_el) {
    res <- mapply(function(x, fmt) sprintf(fmt, x), x = row_el, fmt = nr_fmt)

    if (stars) {
      r <- paste(c(res, sprintf(star_fmt, star(row_el[length(row_el)]))),
        collapse = clp
      )
    } else {
      r <- paste(res, collapse = clp)
    }

    if (tex) r <- paste(r, " \\\\")
    return(r)
  })

  if (tex) names(txt) <- gsub("->", " \\\\textgreater ", names(txt))

  paste(sprintf(txt_fmt, names(txt)), txt, sep = clp)
}

#' Prints a summary of a transfer-entropy result
#'
#' @param object a transfer_entropy
#' @param digits the number of digits to display, defaults to 4
#' @param probs numeric vector of quantiles for the bootstraps
#' @param ... additional arguments, passed to \code{\link[stats]{printCoefmat}}
#'
#' @return invisible the object
#' @export
#'
#' @examples
#' # construct two time-series
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
#' # Calculate Shannon's Transfer Entropy
#' te_result <- transfer_entropy(x, y, nboot = 100)
#'
#' summary(te_result)
summary.transfer_entropy <- function(object, digits = 4,
                                     probs = c(0, 0.25, 0.5, 0.75, 1), ...) {
  cat(sprintf("%s's Transfer Entropy\n\n", fupper(object$entropy)))
  cat("Coefficients:\n")
  printCoefmat(object$coef, ...)

  if (!is.matrix(object$boot)) {
    boot_res <- c(NULL)
  } else {
    quants <- t(apply(object$boot, 1, function(b) quantile(b, probs = probs)))
    rownames(quants) <- c("X->Y", "Y->X")

    probs_nam <- paste0(probs * 100, "%")
    boot_hd_nam <- c("Direction", probs_nam)
    boot_hd_len <- c(
      9,
      rep(
        max(nchar(probs_nam), digits + 2),
        length(probs)
      )
    )

    boot_hd <- paste(mapply(function(l, t) sprintf(sprintf("%%%ss", l), t),
      l = boot_hd_len, t = boot_hd_nam
    ),
    collapse = "  "
    )

    boot_res <- c(
      sprintf(
        "\nBootstrapped TE Quantiles (%s replications):",
        ncol(object$boot)
      ),
      boot_hd,
      textify_mat(quants, digits = digits, width = boot_hd_len, stars = FALSE)
    )
    boot_res <- paste(boot_res, collapse = "\n")
    cat(boot_res, "\n")
  }

  cat(sprintf("\nNumber of Observations: %i%s",
              object$nobs,
              ifelse(object$entropy == "renyi",
                     sprintf("\nQ: %s", object$q),
                     "")))

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
          ifelse(x < 0.1, ".", "")
        )
      )
    )
  )
}


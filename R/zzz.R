#' @useDynLib RTransferEntropy
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile sd printCoefmat
#' @importFrom future plan multisession multicore sequential
NULL

.onAttach <- function(...) {
  set_quiet(FALSE)
}

# first to upper
fupper <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

# returns an s if n > 1 (i.e., sprintf("we have n = %s sample%s", n, mult_s(n)))
mult_s <- function(n) ifelse(n > 1, "s", "")

# wrapper for calc_te and cal_ete that calculates the values
calc_te_ete <- function(restype = "te",
                        x, y, lx = 1, ly = 1, q = 0.1,
                        entropy = "Shannon",
                        shuffles = 100,
                        type = "quantiles",
                        quantiles = c(5, 95),
                        bins = NULL,
                        limits = NULL,
                        burn = 50,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  restype <- tolower(restype)
  if (!restype %in% c("te", "ete")) {
    stop("Internal Error, restype has to be te or ete")
  }

  # Check for unequal length of time series
  if (length(x) != length(y)) {
    stop("x and y must be of same length.")
  }

  # Check that type is specified correctly
  if (!type %in% c("quantiles", "bins", "limits", "q", "b", "l")) {
    stop("type must be either 'quantiles', 'bins' or 'limits'.")
  }

  if (nchar(type) == 1) {
    if (type == "q") {
      type <- "quantiles"
    } else if (type == "b") {
      type <- "bins"
    } else {
      type <- "limits"
    }
  }

  # Check/Restrict number of classes and Markov order/lags
  if (length(quantiles) > 20 || length(bins) > 20 || length(limits) > 20) {
    stop(paste(
      "Number of classes should not exceed 20. Do not expect sensical results",
      "when using too many classes and/or lags."
    ))
  }

  if (lx > 20 || ly > 20) {
    stop(paste(
      "Markov order/number of lags should not exceed 20.",
      "Do not expect sensical results when using too many classes and/or lags."
    ))
  }

  if (lx != ly) {
    warning(paste(
      "Markov order/number of lags should be identical for both time series to",
      "facilitate interpretation of results. Consider setting lx = ly."
    ))
  }

  # Check that transfer entropy measure is specified correctly
  entropy <- tolower(entropy)
  # Allow for specifying the first character only
  if (nchar(entropy) == 1 && entropy %in% c("s", "r")) {
    entropy <- if (entropy == "s") "shannon" else "renyi"
  }

  if (!entropy %in% c("shannon", "renyi")) {
    stop("entropy must be either 'Shannon' or 'Renyi'.")
  }

  # Check that q is between 0 and 1
  if (entropy == "renyi") {
    if (q < 0) {
      stop("q must follow 0 < q < 1")
    } else if (q >= 1) {
      warning(paste(
        "As q-->1, Renyi transfer entropy converges to Shannon transfer",
        "entropy. Shannon transfer entropy is calculated."
      ))
      entropy <- "shannon"
    }
  }

  # Check quantiles
  if (type == "quantiles" && (min(quantiles) < 0 || max(quantiles) > 100)) {
    stop("Quantiles must be between 0 and 100")
  }

  if (type == "quantiles" && max(quantiles) <= 1) {
    warning(paste(
      "Expected quantiles between 0 and 100 but found between 0 and 1,",
      "multiplying by 100."
    ))
    quantiles <- quantiles * 100
  }

  # Remove missing values
  mis_values <- is.na(x) | is.na(y)
  x <- x[!mis_values]
  y <- y[!mis_values]

  if (length(x) == 0) stop("x and y must have non-missing values.")

  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # only calculate the X->Y direction
  if (entropy == "shannon") {
    te <- calc_te_shannon(y, lx = ly, x, ly = lx)
    if (restype == "ete") {
      consty <- shuffle_shannon(
        x = y,
        lx = ly,
        y = x,
        ly = lx,
        shuffles = shuffles
      )
      ete <- te - consty
      ete <- max(0, ete)
    }
    te <- max(0, te)
  } else {
    # RENYI
    te <- calc_te_renyi(y, lx = ly, x, ly = lx, q = q)
    if (restype == "ete") {
      consty <- shuffle_renyi(
        x = y,
        lx = ly,
        y = x,
        ly = lx,
        shuffles = shuffles,
        q = q
      )
      ete <- te - consty
    }
  }

  if (restype == "ete") {
    return(ete)
  } else {
    return(te)
  }
}

#' Daily stock data for 10 stocks from 2000-2017
#'
#' A dataset containing the daily stock returns for 10 stocks and the S&P 500
#' market returns for the time-period 2000-01-04 until 2017-12-29
#'
#' @format A data frame (or data.table if loaded) with 46940 rows and 4 variables:
#' \describe{
#'   \item{date}{date of the observation}
#'   \item{ticker}{ticker of the stock}
#'   \item{ret}{Return of the stock}
#'   \item{sp500}{Return of the S&P 500 stock market index}
#' }
#' @source yahoo finance using \code{\link[quantmod]{getSymbols}}
"stocks"

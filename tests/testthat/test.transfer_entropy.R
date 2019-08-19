# construct two time-series
set.seed(1234567890)
n <- 1000
x <- rep(0, n + 1)
y <- rep(0, n + 1)
for (i in seq(n)) {
  x[i + 1] <- 0.2 * x[i] + rnorm(1, 0, 2)
  y[i + 1] <- x[i] + rnorm(1, 0, 2)
}
x <- x[-1]
y <- y[-1]

#################################
# Shannon Entropy
#################################

context("Shannon's Entropy")

test_that("transfer_entropy shannon is correctly specified", {
  # calc_te X->Y
  res <- calc_te(x, y)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.1120278, tolerance = 1e-6)

  # calc_te Y->X
  res <- calc_te(y, x)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.007642886, tolerance = 1e-6)

  # calc_ete X->Y
  res <- calc_ete(x, y, seed = 1234567890)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.1052938, tolerance = 1e-6)

  # calc_ete Y->X
  res <- calc_ete(y, x, seed = 1234567890)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.001819, tolerance = 1e-6)

  # transfer_entropy
  suppressWarnings({
    res <- transfer_entropy(x, y,
      lx = 1, ly = 1, nboot = 10, quiet = T,
      seed = 12345667
    )
  })

  # check types
  expect_true(is.transfer_entropy(res))

  # we have the all observations saved properly
  expect_equal(res$entropy, "shannon")
  expect_equal(res$obs$x, x)
  expect_equal(res$obs$y, y)
  expect_equal(res$nobs, n)

  # check the coefficients
  coefs <- coef(res)
  expect_true(is.matrix(coefs))
  expect_equal(dim(coefs), c(2, 4))

  # check bootstrapped results
  boot <- res$boot
  expect_true(is.matrix(boot))
  expect_equal(dim(boot), c(2, 10))

  # check values
  exp_coefs <- matrix(
    c(0.112027, 0.007642, 0.104721, 0.002144, 0.003884, 0.002741, 0, 0.3),
    nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})


test_that("transfer_entropy handles missing values", {
  x[10] <- NA
  te <- calc_te(x, y, na.rm = TRUE)
  te2 <- calc_te(x, y, na.rm = FALSE)

  expect_false(is.na(te))
  expect_true(is.na(te2))
})

#################################
# Renyi Entropy
#################################

context("Renyi's Entropy")

test_that("transfer_entropy renyi is correctly specified", {
  # calc_te X->Y
  res <- calc_te(x, y, entropy = "renyi")
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.2530839, tolerance = 1e-6)

  # calc_te Y->X
  res <- calc_te(y, x, entropy = "renyi")
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, 0.02494136, tolerance = 1e-6)

  # calc_ete X->Y
  res <- calc_ete(x, y, seed = 1234567890, entropy = "renyi")
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, -0.072464, tolerance = 1e-6)

  # calc_ete Y->X
  res <- calc_ete(y, x, entropy = "renyi")
  expect_type(res, "double")
  expect_length(res, 1)
  expect_equal(res, -0.169475, tolerance = 1e-6)

  # transfer_entropy
  suppressWarnings({
    res <- transfer_entropy(x, y,
      lx = 1, ly = 1, entropy = "renyi", q = 0.5,
      nboot = 10, quiet = T, seed = 12345667
    )
  })

  # check types
  expect_true(is.transfer_entropy(res))

  # we have the all observations saved properly
  expect_equal(res$entropy, "renyi")
  expect_equal(res$obs$x, x)
  expect_equal(res$obs$y, y)
  expect_equal(res$nobs, n)

  # check the coefficients
  coefs <- coef(res)
  expect_true(is.matrix(coefs))
  expect_equal(dim(coefs), c(2, 4))

  # check bootstrapped results
  boot <- res$boot
  expect_true(is.matrix(boot))
  expect_equal(dim(boot), c(2, 10))

  # check values
  exp_coefs <- matrix(
    c(0.121448, 0.01247, 0.045051, -0.036961, 0.028008, 0.042837, 0, 0.8),
    nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})

#################################
# zoo and xts compatability
#################################

context("zoo & xts compatability")

test_that("Check that transfer_entropy takes zoos and xts", {
  x <- x[1:200]
  y <- y[1:200]
  x.date <- seq(
    from = as.Date("2010-01-01"),
    by = "day",
    length.out = length(x)
  )

  suppressWarnings({
    te_raw <- transfer_entropy(x, y, seed = 123, nboot = 10, quiet = T)
  })

  # ts
  x_ts <- ts(x, start = min(x.date), end = max(x.date))
  y_ts <- ts(y, start = min(x.date), end = max(x.date))
  suppressWarnings({
    te_ts <- transfer_entropy(x_ts, y_ts, seed = 123, nboot = 10, quiet = T)
  })
  expect_equal(te_raw, te_ts)

  # zoo
  x_zoo <- zoo::zoo(x, x.date)
  y_zoo <- zoo::zoo(y, x.date)

  suppressWarnings({
    te_zoo <- transfer_entropy(x_zoo, y_zoo, seed = 123, nboot = 10, quiet = T)
  })
  expect_equal(te_raw, te_zoo)

  # xts
  x_xts <- xts::xts(x, x.date)
  y_xts <- xts::xts(y, x.date)
  suppressWarnings({
    te_xts <- transfer_entropy(x_xts, y_xts, seed = 123, nboot = 10, quiet = T)
  })
  expect_equal(te_raw, te_xts)
})

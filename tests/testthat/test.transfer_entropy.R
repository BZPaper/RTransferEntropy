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


test_that("te_result is correctly specified", {
  context("Shannon")

  suppressWarnings({
    res <- transfer_entropy(x, y, lx = 1, ly = 1, nboot = 10, quiet = T,
                            seed = 12345667)
  })

  context("check types")
  expect_true(is.TEResult(res))

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

  context("check values")
  exp_coefs <- matrix(
    c(0.0157868988194173, 0.00173401913231719, 0.0126125335635818,
      0, 0.00197686371110982, 0.00130962987420229, 0, 0.9), nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})


#################################
# Renyi Entropy
#################################

test_that("te_result is correctly specified", {
  context("Renyi")

  suppressWarnings({
    res <- transfer_entropy(x, y, lx = 1, ly = 1, entropy = "renyi", q = 0.5,
                            nboot = 10, quiet = T, seed = 12345667)
  })

  context("check types")
  expect_true(is.TEResult(res))

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

  context("check values")
  exp_coefs <- matrix(
    c(0.0607475858536458, 0.0112538569883853, 0.0120269808040249,
      0, 0.021012313171769, 0.0191147381370581, 0.2, 1),
    nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})

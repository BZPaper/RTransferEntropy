# construct two time-series
set.seed(1234567890)
n <- 10000
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

  suppressWarnings(
    res <- transfer_entropy(x, y, lx = 1, ly = 1, nboot = 10, quiet = T,
                            seed = 12345667, cl = 1)
  )

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
    c(0.0996107279340461, 0.00128027810811011, 0.0986370632185543,
      0.000392488584173327, 0.000386084463301625, 0.000489447148659146,
      0, 0.2), nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})


#################################
# Renyi Entropy
#################################

test_that("te_result is correctly specified", {
  context("Renyi")

  suppressWarnings(
    res <- transfer_entropy(x, y, lx = 1, ly = 1, entropy = "renyi", q = 0.5,
                            nboot = 10, quiet = T, seed = 12345667, cl = 1)
  )
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
    c(0.0941334355990759, 0.0298813788644767, 0.0765819430684522,
      0.0149377788357564, 0.01477246415148, 0.0098556260764617, 0,
      0.1),
    nrow = 2, ncol = 4,
    dimnames = list(c("X->Y", "Y->X"), c("te", "ete", "se", "p-value"))
  )
  expect_equal(coefs, exp_coefs, tolerance = 1e-6)
})

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
  res <- transfer_entropy(x, y, lx = 1, ly = 1, nboot = 10, quiet = T,
                          seed = 12345667)

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
  exp_coefs <- matrix(c(0.0996107279340461, 0.00128027810811011,
                        0.0987866909219102, 0.000472164215081831,
                        0.000108401564915032, 0.000119239167199667,
                        0, 0), nrow = 2, ncol = 4,
                      dimnames = list(c("X->Y", "Y->X"),
                                      c("te", "ete", "se", "p-value")))
  expect_equal(coefs, exp_coefs)
})


#################################
# Renyi Entropy
#################################

test_that("te_result is correctly specified", {
  context("Renyi")
  res <- transfer_entropy(x, y, lx = 1, ly = 1, entropy = "renyi", q = 0.5,
                          nboot = 10, quiet = T, seed = 12345667)

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
  exp_coefs <- matrix(c(0.0941334355990759, 0.0298813788644767,
                        0.0808040232549858, 0.0136590056511369,
                        0.00401946422132985, 0.00316381604679907,
                        0, 0), nrow = 2, ncol = 4,
                      dimnames = list(c("X->Y", "Y->X"),
                                      c("te", "ete", "se", "p-value")))
  expect_equal(coefs, exp_coefs)
})

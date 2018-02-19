library(RTransferEntropy)
library(testthat)

# tests if x is numeric (allowing for NAs) and of length 1
expect_numeric_1 <- function(x) (is.numeric(x) | is.na(x)) && length(x) == 1

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

res <- transfer_entropy(x, lx = 1, y, ly = 1, quiet = T, seed = 12345667)

context("Result Type")
test_that("te_result is correctly specified", {
  expect_true(is.TEResult(res))

  # we have the all observations saved properly
  expect_equal(res$entropy, "shannon")
  expect_equal(res$obs$x, x)
  expect_equal(res$obs$y, y)
  expect_equal(res$nobs, n)

  # we have only one result for each te, ete, sete, pete
  expect_true(expect_numeric_1(res$te_xy))
  expect_true(expect_numeric_1(res$te_yx))
  expect_true(expect_numeric_1(res$ete_xy))
  expect_true(expect_numeric_1(res$ete_yx))
  expect_true(expect_numeric_1(res$se_xy))
  expect_true(expect_numeric_1(res$se_yx))
  expect_true(expect_numeric_1(res$p_xy))
  expect_true(expect_numeric_1(res$p_yx))
})

context("Result Value")
test_that("te_result have certain values for the seed 12345667", {
  expect_equal(res$te_xy, 0.09961073)
  expect_equal(res$te_yx, 0.001280278)
  expect_equal(res$ete_xy, 0.09872741)
  expect_equal(res$ete_yx, 0.0003214504)
  expect_equal(res$se_xy, 0.0001084016)
  expect_equal(res$se_yx, 0.0001192392)
  expect_equal(res$p_xy, 0)
  expect_equal(res$p_yx, 1)
})


#################################
# Renyi Entropy
#################################

res <- transfer_entropy(x, lx = 1, y, ly = 1, entropy = "renyi", q = 0.5, cl = 1, quiet = T, seed = 12345667)


context("Result Type")
test_that("te_result is correctly specified", {
  expect_true(is.TEResult(res))

  # we have the all observations saved properly
  expect_equal(res$entropy, "renyi")
  expect_equal(res$obs$x, x)
  expect_equal(res$obs$y, y)
  expect_equal(res$nobs, n)

  # we have only one result for each te, ete, sete, pete
  expect_true(expect_numeric_1(res$te_xy))
  expect_true(expect_numeric_1(res$te_yx))
  expect_true(expect_numeric_1(res$ete_xy))
  expect_true(expect_numeric_1(res$ete_yx))
  expect_true(expect_numeric_1(res$se_xy))
  expect_true(expect_numeric_1(res$se_yx))
  expect_true(expect_numeric_1(res$p_xy))
  expect_true(expect_numeric_1(res$p_yx))
})

context("Result Value")
test_that("te_result have certain values for the seed 12345667", {
  expect_equal(res$te_xy, 0.09413344)
  expect_equal(res$te_yx, 0.02988138)
  expect_equal(res$ete_xy, 0.07991739)
  expect_equal(res$ete_yx, 0.0104142)
  expect_equal(res$se_xy, 0.003887758)
  expect_equal(res$se_yx, 0.002233663)
  expect_equal(res$p_xy, 0)
  expect_equal(res$p_yx, 1)
})

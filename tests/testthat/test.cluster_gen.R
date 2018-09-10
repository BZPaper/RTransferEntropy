x <- c(1, 2, 2, 1, 1, 2)

context("Check cluster_gen")

test_that("cluster_gen with no lag", {
  cls <- RTransferEntropy:::cluster_gen(x, prog = FALSE)
  frqs_exp <- c("1" = 0.6, "2" = 0.4)
  expect_equal(cls, frqs_exp)

  cls <- RTransferEntropy:::cluster_gen(x, prog = TRUE)
  frqs_exp <- c("1 1" = 0.2, "1 2" = 0.4, "2 1" = 0.2, "2 2" = 0.2)
  expect_equal(cls, frqs_exp)
})

test_that("cluster_gen with lag of x of 2", {
  cls <- RTransferEntropy:::cluster_gen(x, lx = 2, prog = FALSE)
  frqs_exp <- c("1 1" = 0.25, "1 2" = 0.25, "2 1" = 0.25, "2 2" = 0.25)
  expect_equal(cls, frqs_exp)

  cls <- RTransferEntropy:::cluster_gen(x, lx = 2, prog = TRUE)
  frqs_exp <- c("1 1 2" = 0.25, "1 2 2" = 0.25, "2 1 1" = 0.25, "2 2 1" = 0.25)
  expect_equal(cls, frqs_exp)
})

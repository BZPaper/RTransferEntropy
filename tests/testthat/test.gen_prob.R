set.seed(123)
x <- sample(1:3, 10, replace = T)

test_that("gen_prob is correctly specified", {

  context("gen_prob")
  gen <- RTransferEntropy:::gen_prob(x, 1)

  context("Check types")
  expect_true(is.list(gen))
  expect_equal(length(gen), 2)
  expect_equal(sum(gen$px), 1)

  expect_true(is.list(gen$transprob))
  expect_equal(length(gen$transprob), 3)
  # in transprob all values sum up to 1
  expect_equal(sapply(gen$transprob, sum), c("1" = 1, "2" = 1, "3" = 1))


  context("Check values")
  # gen$px is a table... not a named vector... as it is created using the table-function
  expect_equal(as.numeric(gen$px), c(0.2, 0.4, 0.4))
  expect_equal(gen$transprob[[1]], c("12" = 0.5, "13" = 0.5))
  expect_equal(gen$transprob[[2]], c("22" = 1/3, "23" = 2/3))
  expect_equal(gen$transprob[[3]], c("31" = 0.25, "32" = 0.5, "33" = 0.25))
})

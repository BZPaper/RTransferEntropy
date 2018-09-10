context("Check calculate_transisition_probabilities")

test_that("Check calculate_transition_probabilities for complete cases", {
  x <- c(1, 2, 3, 2, 1)

  res <- RTransferEntropy:::calculate_transition_probabilities(x, 1)

  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_equal(names(res), c("1", "2", "3"))
  expect_equal(res[[1]], c("1 2" = 1))
  expect_equal(res[[2]], c("2 1" = 0.5, "2 3" = 0.5))
  expect_equal(res[[3]], c("3 2" = 1))
})

test_that("Check calculate_transition_probabilities incomplete cases", {
  x <- c(1, 1, 3, 3, 3)

  res <- RTransferEntropy:::calculate_transition_probabilities(x, 1)

  expect_equal(length(res), 3)
  expect_equal(names(res), c("1", "2", "3"))
  expect_equal(res[[1]], c("1 1" = 0.5, "1 3" = 0.5))
  expect_equal(res[[2]], structure(numeric(0), .Names = character(0)))
  expect_equal(res[[3]], c("3 3" = 1))
})

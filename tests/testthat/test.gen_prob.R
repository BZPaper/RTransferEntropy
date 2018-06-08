
test_that("Check gen_prob for lx = 1", {
  set.seed(123)
  x <- sample(1:3, 10, replace = T)

  context("Check lag = 1")
  gen <- RTransferEntropy:::gen_prob(x, 1)

  context("Check types")
  expect_true(is.list(gen))
  expect_equal(length(gen), 2)
  expect_equal(sum(gen$px), 1)

  expect_true(is.list(gen$transprob))
  expect_equal(length(gen$transprob), 3)
  # in transprob all values sum up to 1
  expect_equal(sapply(gen$transprob, sum), c("1" = 1, "2" = 1, "3" = 1))

  context("Check Transition Probabilities")
  expect_equal(as.numeric(gen$px), c(0.2, 0.4, 0.4))
  expect_equal(gen$transprob[[1]], c("12" = 0.5, "13" = 0.5))
  expect_equal(gen$transprob[[2]], c("22" = 1 / 3, "23" = 2 / 3))
  expect_equal(gen$transprob[[3]], c("31" = 0.25, "32" = 0.5, "33" = 0.25))
})


test_that("Check gen_prob for lx = 1s", {
  set.seed(123)
  k <- 4
  n <- 1000
  x <- sample(1:k, n, replace = T)

  context("Check lag > 1")
  gen <- RTransferEntropy:::gen_prob(x, 2)

  expect_true(is.list(gen))
  expect_equal(length(gen), 2)
  expect_equal(sum(gen$px), 1)

  expect_true(is.list(gen$transprob))
  expect_equal(length(gen$transprob), k)
  # in transprob all values sum up to 1
  exp_sum <- rep(1, k)
  names(exp_sum) <- 1:k
  expect_equal(sapply(gen$transprob, sum), exp_sum)

  context("Check Transition Probabilities")
  expect_equal(as.numeric(gen$px), c(0.247, 0.260, 0.245, 0.248))
  prob1 <- c(0.044715, 0.069106, 0.056911, 0.081301, 0.056911, 0.089431,
             0.060976, 0.069106, 0.052846, 0.081301, 0.036585, 0.056911,
             0.044715, 0.077236, 0.069106, 0.052846)
  names(prob1) <- c("111", "112", "113", "114", "121", "122", "123", "124",
                    "131", "132", "133", "134", "141", "142", "143", "144")
  expect_equal(gen$transprob[[1]], prob1, tolerance = 1e-6)

  prob2 <- c(0.080769, 0.026923, 0.080769, 0.05, 0.080769, 0.065385, 0.065385,
             0.057692, 0.073077, 0.038462, 0.057692, 0.069231, 0.080769,
             0.042308, 0.061538, 0.069231)
  names(prob2) <- c("211", "212", "213", "214", "221", "222", "223", "224",
                    "231", "232", "233", "234", "241", "242", "243", "244")
  expect_equal(gen$transprob[[2]], prob2, tolerance = 1e-6)

  prob3 <- c(0.04918, 0.077869, 0.040984, 0.057377, 0.053279, 0.04918, 0.07377,
             0.086066, 0.061475, 0.065574, 0.069672, 0.053279, 0.065574,
             0.053279, 0.086066, 0.057377)
  names(prob3) <- c("311", "312", "313", "314", "321", "322", "323", "324",
                    "331", "332", "333", "334", "341", "342", "343", "344")
  expect_equal(gen$transprob[[3]], prob3, tolerance = 1e-6)


  prob4 <- c(0.072581, 0.100806, 0.044355, 0.052419, 0.056452, 0.076613,
             0.048387, 0.048387, 0.03629, 0.072581, 0.080645, 0.076613,
             0.076613, 0.056452, 0.048387, 0.052419)
  names(prob4) <- c("411", "412", "413", "414", "421", "422", "423", "424",
                    "431", "432", "433", "434", "441", "442", "443", "444")
  expect_equal(gen$transprob[[4]], prob4, tolerance = 1e-6)
})

context("classes.R")

test_that("rsi works", {
  expect_true(as.rsi("S") < as.rsi("I"))
  expect_true(as.rsi("I") < as.rsi("R"))
  expect_true(as.rsi("R") > as.rsi("S"))
  expect_true(is.rsi(as.rsi("S")))

  # print plots, should not raise errors
  barplot(as.rsi(c("S", "I", "R")))
  plot(as.rsi(c("S", "I", "R")))
  print(as.rsi(c("S", "I", "R")))

  expect_equal(suppressWarnings(as.logical(as.rsi("INVALID VALUE"))), NA)

  expect_equal(summary(as.rsi(c("S", "R"))), c("Mode" = 'rsi',
                                               "<NA>" = "0",
                                               "Sum S" = "1",
                                               "Sum IR" = "1",
                                               "Sum R" = "1",
                                               "Sum I" = "0"))
})

test_that("mic works", {
  expect_true(as.mic(8) == as.mic("8"))
  expect_true(as.mic("1") > as.mic("<=0.0625"))
  expect_true(as.mic("1") < as.mic(">=32"))
  expect_true(is.mic(as.mic(8)))

  expect_equal(as.double(as.mic(">=32")), 32)
  expect_equal(as.integer(as.mic(">=32")), 32)
  expect_equal(suppressWarnings(as.logical(as.mic("INVALID VALUE"))), NA)

  # print plots, should not raise errors
  barplot(as.mic(c(1, 2, 4, 8)))
  plot(as.mic(c(1, 2, 4, 8)))
  print(as.mic(c(1, 2, 4, 8)))

  expect_equal(summary(as.mic(c(2, 8))), c("Mode" = 'mic',
                                           "<NA>" = "0",
                                           "Min." = "2",
                                           "Max." = "8"))
})

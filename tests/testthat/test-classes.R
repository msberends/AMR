context("classes.R")

test_that("rsi works", {
  expect_true(as.rsi("S") < as.rsi("I"))
  expect_true(as.rsi("I") < as.rsi("R"))
  expect_true(as.rsi("R") > as.rsi("S"))
  expect_true(is.rsi(as.rsi("S")))
  
  expect_equal(suppressWarnings(as.logical(as.rsi("INVALID VALUE"))), NA)
  
  expect_equal(class(barplot(as.rsi(c("S", "I", "R")))), "numeric")
  
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
  
  expect_equal(class(plot(as.mic(c(1, 2, 4, 8)))), "numeric")
  
  expect_equal(summary(as.mic(c(2, 8))), c("Mode" = 'mic',
                                           "<NA>" = "0",
                                           "Min." = "2",
                                           "Max." = "8"))
})

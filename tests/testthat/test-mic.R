context("mic.R")

test_that("mic works", {
  expect_true(as.mic(8) == as.mic("8"))
  expect_true(as.mic("1") > as.mic("<=0.0625"))
  expect_true(as.mic("1") < as.mic(">=32"))
  expect_true(is.mic(as.mic(8)))

  expect_equal(as.double(as.mic(">=32")), 32)
  expect_equal(as.numeric(as.mic(">=32")), 32)
  expect_equal(as.integer(as.mic(">=32")), 32)
  expect_equal(suppressWarnings(as.logical(as.mic("INVALID VALUE"))), NA)

  # all levels should be valid MICs
  expect_silent(as.mic(levels(as.mic(1))))

  expect_warning(as.mic("INVALID VALUE"))

  # print plots, should not raise errors
  barplot(as.mic(c(1, 2, 4, 8)))
  plot(as.mic(c(1, 2, 4, 8)))
  print(as.mic(c(1, 2, 4, 8)))

  expect_equal(summary(as.mic(c(2, 8))), c("Mode" = 'mic',
                                           "<NA>" = "0",
                                           "Min." = "2",
                                           "Max." = "8"))
})

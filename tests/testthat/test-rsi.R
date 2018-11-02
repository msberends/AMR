context("rsi.R")

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
                                               "-Sum R" = "1",
                                               "-Sum I" = "0"))

  expect_identical(as.logical(lapply(septic_patients, is.rsi.eligible)),
                   rep(FALSE, length(septic_patients)))

  library(dplyr)
  # 40 rsi columns
  expect_equal(septic_patients %>%
                 mutate_at(vars(peni:rifa), as.character) %>%
                 lapply(is.rsi.eligible) %>%
                 as.logical() %>%
                 sum(),
               40)

})

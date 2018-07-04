context("p.symbol.R")

test_that("P symbol works", {
  expect_identical(p.symbol(c(0.001, 0.01, 0.05, 0.1, 1, NA, 3)),
                   c("***", "**", "*", ".", " ", NA, NA))
})

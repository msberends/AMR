context("misc.R")

test_that("`like` works", {
  expect_true("test" %like% "^t")
  expect_true("test" %like% "test")
  expect_true("test" %like% "TEST")
  expect_true(as.factor("test") %like% "TEST")
})

test_that("percentages works", {
  expect_equal(percent(0.25), "25%")
  expect_equal(percent(0.5), "50%")
  expect_equal(percent(0.500, force_zero = TRUE), "50.0%")
  expect_equal(percent(0.1234), "12.3%")
})


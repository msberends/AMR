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

test_that("size format works", {
  expect_equal(size_humanreadable(123456), "121 kB")
})

test_that("functions missing in older R versions work", {
  expect_equal(strrep2("A", 5), "AAAAA")
  expect_equal(strrep2(c("A", "B"), c(5, 2)), c("AAAAA", "BB"))
  expect_equal(trimws(" test "), "test")
  expect_equal(trimws(" test ", "l"), "test ")
  expect_equal(trimws(" test ", "r"), " test")
})

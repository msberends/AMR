context("misc.R")

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
  expect_equal(strrep("A", 5), "AAAAA")
  expect_equal(strrep(c("A", "B"), c(5, 2)), c("AAAAA", "BB"))
  expect_equal(trimws(" test "), "test")
  expect_equal(trimws(" test ", "l"), "test ")
  expect_equal(trimws(" test ", "r"), " test")
})

test_that("generic dates work", {
  expect_equal(date_generic("yyyy-mm-dd"), "%Y-%m-%d")
  expect_equal(date_generic("dddd d mmmm yyyy"), "%A %e %B %Y")
})

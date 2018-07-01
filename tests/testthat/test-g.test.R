context("g.test.R")

test_that("G-test works", {

  # example 1: clearfield rice vs. red rice
  x <- c(772, 1611, 737)
  x.expected <- vector2ratio(x, ratio = "1:2:1")
  expect_equal(g.test(x, x.expected),
               expected = 0.12574,
               tolerance = 0.00001)

  # example 2: red crossbills
  x <- c(1752, 1895)
  x.expected <- vector2ratio(x, ratio = c(1, 1))
  expect_equal(g.test(x, x.expected),
               expected = 0.01787343,
               tolerance = 0.00000001)

})

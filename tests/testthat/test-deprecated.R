context("deprecated.R")

test_that("deprecated functions work", {

  expect_identical(is.mo(as.mo("esco")), suppressWarnings(is.bactid(as.bactid("esco"))))
  expect_warning(identical(is.mo(as.mo("esco")), is.bactid(as.bactid("esco"))))

  expect_identical(as.mo("esco"), suppressWarnings(guess_bactid("esco")))

  expect_error(suppressWarnings(ratio("A")))
  expect_error(suppressWarnings(ratio(1, ratio = "abc")))
  expect_error(suppressWarnings(ratio(c(1, 2), ratio = c(1, 2, 3))))
  expect_warning(ratio(c(772, 1611, 737), ratio = "1:2:1"))
  expect_identical(suppressWarnings(ratio(c(772, 1611, 737), ratio = "1:2:1")), c(780, 1560,  780))
  expect_identical(suppressWarnings(ratio(c(1752, 1895), ratio = c(1, 1))), c(1823.5, 1823.5))

  old_mo <- "ESCCOL"
  class(old_mo) <- "bactid"
  # print
  expect_output(print(old_mo))
  # test data.frame and pull
  expect_equal(as.character(dplyr::pull(data.frame(test = old_mo), test)), "ESCCOL")

})

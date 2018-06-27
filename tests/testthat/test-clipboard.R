context("clipboard.R")

test_that("clipboard works", {
  skip_if_not(clipr::clipr_available())

  clipboard_export(antibiotics)
  expect_identical(antibiotics,
                   clipboard_import())
})

context("clipboard.R")

test_that("clipboard works", {
  skip_if_not(clipr::clipr_available())

  clipboard_export(antibiotics)
  expect_identical(antibiotics,
                   clipboard_import(date_format = "yyyy-mm-dd"))

  clipboard_export(septic_patients[1:100,])
  expect_identical(tbl_parse_guess(septic_patients[1:100,]),
                   clipboard_import(guess_col_types = TRUE))
})

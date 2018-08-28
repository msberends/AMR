context("clipboard.R")

test_that("clipboard works", {

  if (!clipr::clipr_available() & Sys.info()['sysname'] == "Linux") {
    # try to support on X11, by setting the R variable DISPLAY
    Sys.setenv(DISPLAY = "localhost:10.0")
  }

  skip_if_not(clipr::clipr_available())

  # clipboard_export(antibiotics)
  # imp <- clipboard_import(guess_col_types = FALSE,
  #                         stringsAsFactors = FALSE)
  # expect_identical(as.data.frame(antibiotics, stringsAsFactors = FALSE),
  #                  imp)

  clipboard_export(septic_patients[1:100,])
  imp <- clipboard_import(guess_col_types = TRUE,
                          stringsAsFactors = FALSE)
  expect_identical(as.data.frame(tbl_parse_guess(septic_patients[1:100,]),
                                 stringsAsFactors = FALSE),
                   imp)
})

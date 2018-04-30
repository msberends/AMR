context("clipboard.R")

test_that("clipboard works", {
  if (grepl(Sys.info()['sysname'], "windows", ignore.case = TRUE)) {
    t1 <<- AMR::antibiotics # why is the <<- needed? Won't work without it...
    clipboard_export(t1, info = FALSE)
    t2 <- clipboard_import()
    expect_equal(t1, t2)
  } else {
    expect_equal(TRUE, TRUE)
  }
})

context("clipboard.R")

test_that("clipboard works", {
  t1 <<- AMR::antibiotics # why is the <<- needed? Won't work without it...
  clipboard_export(t1, info = FALSE)
  t2 <- clipboard_import()
  skip_if(is.null(t1) | is.null(t2), message = "No clipboard content found: skipping.")
  expect_equal(t1, t2)
})

context("clipboard.R")

test_that("clipboard works", {
  # why is the <<- needed? Won't work without it...
  t1 <<- AMR::antibiotics
  clipboard_export(t1, info = FALSE)
  t2 <- clipboard_import()
  expect_equal(t1, t2)
})

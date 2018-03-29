context("clipboard.R")

test_that("clipboard works", {
  skip_on_os(c("linux", "solaris"))
  t1 <<- AMR::antibiotics # why is the <<- needed? Won't work without it...
  clipboard_export(t1, info = FALSE)
  t2 <- clipboard_import()
  expect_equal(t1, t2)
})

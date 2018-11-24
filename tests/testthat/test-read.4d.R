context("read.4d.R")

test_that("read 4D works", {

  library(dplyr)
  test1 <- data.frame(Patientnr = "ABC",
                      MV = "M",
                      Monsternr = "0123",
                      Afnamedat = "10-11-12",
                      Bepaling = "bk",
                      Afd. = "ABC",
                      Spec = "ABC",
                      Matbijz. = "ABC",
                      Mat = "ABC",
                      Mocode = "esccol",
                      PENI = "R",
                      stringsAsFactors = FALSE)
  tf <- tempfile()
  write.table(test1, file = tf, quote = F, sep = "\t")

  x <- read.4D(tf, skip = 0, info = TRUE)
  unlink(tf)

  expect_equal(ncol(x), 11)
  expect_equal(class(x$date_received), "Date")
  expect_equal(class(x$mo), "mo")
  expect_equal(as.character(x$mo), "B_ESCHR_COL")
  expect_equal(is.rsi(x$peni), TRUE)

})

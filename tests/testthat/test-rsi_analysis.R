context("rsi_analysis.R")

test_that("rsi works", {
  # amox resistance in `septic_patients` should be around 53.86%
  amox_R <- septic_patients %>% summarise(amox = rsi(amox)) %>% pull(amox)
  expect_equal(amox_R, 0.5386, tolerance = 0.0001)
  expect_equal(rsi_df(septic_patients, 
                      ab = "amox",
                      info = FALSE), 0.5386, tolerance = 0.0001)
  # and pita+genta susceptibility around 98.09%
  expect_equal(rsi_df(septic_patients,
                      ab = c("pita", "gent"), 
                      interpretation = "S", 
                      info = FALSE), 0.9809, tolerance = 0.0001)
})

test_that("prediction of rsi works", {
  amox_R <- rsi_predict(tbl = septic_patients[which(septic_patients$bactid == "ESCCOL"),], 
                        col_ab = "amox",
                        col_date = "date",
                        info = FALSE)
  amox_R <- amox_R %>% pull("probR")
  # amox resistance will decrease according to `septic_patients`
  expect_true(amox_R[2] > amox_R[20])
  expect_error(rsi_predict(tbl = septic_patients[which(septic_patients$bactid == "ESCCOL"),], 
                           model = "INVALID MODEL",
                           col_ab = "amox",
                           col_date = "date",
                           info = FALSE))
  expect_error(rsi_predict(tbl = septic_patients[which(septic_patients$bactid == "ESCCOL"),], 
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = FALSE))
  expect_error(rsi_predict(tbl = septic_patients[which(septic_patients$bactid == "ESCCOL"),], 
                           col_ab = "amox",
                           col_date = "NOT EXISTING COLUMN",
                           info = FALSE))
})

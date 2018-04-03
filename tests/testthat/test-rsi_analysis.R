context("rsi_analysis.R")

test_that("rsi works", {
  # amox resistance in `septic_patients` should be around 53.86%
  expect_equal(rsi(septic_patients$amox), 0.5386, tolerance = 0.0001)
  expect_equal(rsi(septic_patients$amox), 0.5386, tolerance = 0.0001)
  expect_equal(rsi_df(septic_patients,
                      ab = "amox",
                      info = FALSE),
               0.5386,
               tolerance = 0.0001)
  # pita+genta susceptibility around 98.09%
  expect_equal(rsi(septic_patients$pita,
                   septic_patients$gent,
                   interpretation = "S",
                   info = TRUE),
               0.9809,
               tolerance = 0.0001)
  expect_equal(rsi_df(septic_patients,
                      ab = c("pita", "gent"),
                      interpretation = "S",
                      info = FALSE),
               0.9809,
               tolerance = 0.0001)
  # mero+pita+genta susceptibility around 98.58%
  expect_equal(rsi_df(septic_patients,
                      ab = c("mero", "pita", "gent"),
                      interpretation = "IS",
                      info = FALSE),
               0.9858,
               tolerance = 0.0001)
})

test_that("prediction of rsi works", {
  amox_R <- septic_patients %>%
    filter(bactid == "ESCCOL") %>%
    rsi_predict(col_ab = "amox",
                col_date = "date",
                info = FALSE) %>%
    pull("probR")
  # amox resistance will decrease using dataset `septic_patients`
  expect_true(amox_R[2] > amox_R[20])

  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "binomial",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "loglin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "lin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))

  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           model = "INVALID MODEL",
                           col_ab = "amox",
                           col_date = "date",
                           info = FALSE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = FALSE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "amox",
                           col_date = "NOT EXISTING COLUMN",
                           info = FALSE))
})

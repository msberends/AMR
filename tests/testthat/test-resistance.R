context("resistance.R")

test_that("resistance works", {
  # amox resistance in `septic_patients` should be around 57.56%
  expect_equal(resistance(septic_patients$amox), 0.5756, tolerance = 0.0001)
  expect_equal(susceptibility(septic_patients$amox), 1 - 0.5756, tolerance = 0.0001)

  # pita+genta susceptibility around 98.09%
  expect_equal(susceptibility(septic_patients$pita,
                              septic_patients$gent),
               0.9809,
               tolerance = 0.0001)
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S")),
               0.9809,
               tolerance = 0.0001)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_p = susceptibility(cipr, as_percent = TRUE),
                           cipro_n = n_rsi(cipr),
                           genta_p = susceptibility(gent, as_percent = TRUE),
                           genta_n = n_rsi(gent),
                           combination_p = susceptibility(cipr, gent, as_percent = TRUE),
                           combination_n = n_rsi(cipr, gent)) %>%
                 pull(combination_n),
               c(138, 474, 170, 464, 183))
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

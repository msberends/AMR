context("rsi_IR.R")

test_that("resistance works", {
  # check shortcuts
  expect_equal(resistance(septic_patients$amox, include_I = TRUE),
               rsi_IR(septic_patients$amox))
  expect_equal(resistance(septic_patients$amox, include_I = FALSE),
               rsi_R(septic_patients$amox))
  expect_equal(intermediate(septic_patients$amox),
               rsi_I(septic_patients$amox))
  expect_equal(susceptibility(septic_patients$amox, include_I = TRUE),
               rsi_SI(septic_patients$amox))
  expect_equal(susceptibility(septic_patients$amox, include_I = FALSE),
               rsi_S(septic_patients$amox))

  # amox resistance in `septic_patients`
  expect_equal(rsi_R(septic_patients$amox), 0.6603, tolerance = 0.0001)
  expect_equal(rsi_I(septic_patients$amox), 0.0030, tolerance = 0.0001)
  expect_equal(1 - rsi_R(septic_patients$amox) - rsi_I(septic_patients$amox),
               rsi_S(septic_patients$amox))

  # pita+genta susceptibility around 98.09%
  expect_equal(susceptibility(septic_patients$pita,
                              septic_patients$gent),
               0.9535,
               tolerance = 0.0001)
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S")),
               0.9535,
               tolerance = 0.0001)

  # percentages
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(R = rsi_R(cipr, as_percent = TRUE),
                           I = rsi_I(cipr, as_percent = TRUE),
                           S = rsi_S(cipr, as_percent = TRUE),
                           n = rsi_n(cipr),
                           total = n()) %>%
                 pull(n) %>%
                 sum(),
               1404)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_p = susceptibility(cipr, as_percent = TRUE),
                           cipro_n = n_rsi(cipr),
                           genta_p = susceptibility(gent, as_percent = TRUE),
                           genta_n = n_rsi(gent),
                           combination_p = susceptibility(cipr, gent, as_percent = TRUE),
                           combination_n = rsi_n(cipr, gent)) %>%
                 pull(combination_n),
               c(202, 482, 201, 499))

  expect_warning(resistance(as.character(septic_patients$amcl)))
  expect_warning(susceptibility(as.character(septic_patients$amcl)))
  expect_warning(susceptibility(as.character(septic_patients$amcl,
                                             septic_patients$gent)))


  # check for errors
  expect_error(rsi_IR(septic_patients %>% select(amox, amcl)))
  expect_error(rsi_IR("test", minimum = "test"))
  expect_error(rsi_IR("test", as_percent = "test"))
  expect_error(rsi_I(septic_patients %>% select(amox, amcl)))
  expect_error(rsi_I("test", minimum = "test"))
  expect_error(rsi_I("test", as_percent = "test"))
  expect_error(rsi_S("test", minimum = "test"))
  expect_error(rsi_S("test", as_percent = "test"))
  expect_error(rsi_S(septic_patients %>% select(amox, amcl)))
  expect_error(rsi_S("R", septic_patients %>% select(amox, amcl)))

  # check too low amount of isolates
  expect_identical(rsi_R(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)
  expect_identical(rsi_I(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)
  expect_identical(rsi_S(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)

  # warning for speed loss
  expect_warning(rsi_R(as.character(septic_patients$gent)))
  expect_warning(rsi_I(as.character(septic_patients$gent)))
  expect_warning(rsi_S(septic_patients$amcl, as.character(septic_patients$gent)))

})

test_that("old rsi works", {
  # amox resistance in `septic_patients` should be around 66.33%
  expect_equal(suppressWarnings(rsi(septic_patients$amox)), 0.6633, tolerance = 0.0001)
  expect_equal(suppressWarnings(rsi(septic_patients$amox, interpretation = "S")), 1 - 0.6633, tolerance = 0.0001)
  expect_equal(rsi_df(septic_patients,
                      ab = "amox",
                      info = TRUE),
               0.6633,
               tolerance = 0.0001)
  # pita+genta susceptibility around 98.09%
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S",
                                    info = TRUE)),
               0.9535,
               tolerance = 0.0001)
  expect_equal(rsi_df(septic_patients,
                      ab = c("pita", "gent"),
                      interpretation = "S",
                      info = TRUE),
               0.9535,
               tolerance = 0.0001)
  # more than 2 not allowed
  expect_error(rsi_df(septic_patients,
                      ab = c("mero", "pita", "gent"),
                      interpretation = "IS",
                      info = TRUE))

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_S = suppressWarnings(rsi(cipr, interpretation = "S",
                                                          as_percent = TRUE, warning = FALSE)),
                           cipro_n = n_rsi(cipr),
                           genta_S = suppressWarnings(rsi(gent, interpretation = "S",
                                                          as_percent = TRUE, warning = FALSE)),
                           genta_n = n_rsi(gent),
                           combination_S = suppressWarnings(rsi(cipr, gent, interpretation = "S",
                                                                as_percent = TRUE, warning = FALSE)),
                           combination_n = n_rsi(cipr, gent)) %>%
                 pull(combination_n),
               c(202, 482, 201, 499))
})

test_that("prediction of rsi works", {
  amox_R <- septic_patients %>%
    filter(bactid == "ESCCOL") %>%
    rsi_predict(col_ab = "amox",
                col_date = "date",
                minimum = 10,
                info = TRUE) %>%
    pull("value")
  # amox resistance will increase according to data set `septic_patients`
  expect_true(amox_R[3] < amox_R[20])

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
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "amox",
                           col_date = "NOT EXISTING COLUMN",
                           info = TRUE))
  # almost all E. coli are mero S in the Netherlands :)
  expect_error(resistance_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                                  col_ab = "mero",
                                  col_date = "date",
                                  info = TRUE))
})

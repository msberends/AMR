# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

context("portion.R")

test_that("portions works", {
  # amox resistance in `septic_patients`
  expect_equal(portion_R(septic_patients$amox), 0.5827645, tolerance = 0.0001)
  expect_equal(portion_I(septic_patients$amox), 0.0025597, tolerance = 0.0001)
  expect_equal(1 - portion_R(septic_patients$amox) - portion_I(septic_patients$amox),
               portion_S(septic_patients$amox))
  expect_equal(portion_R(septic_patients$amox) + portion_I(septic_patients$amox),
               portion_IR(septic_patients$amox))
  expect_equal(portion_S(septic_patients$amox) + portion_I(septic_patients$amox),
               portion_SI(septic_patients$amox))

  expect_equal(septic_patients %>% portion_S(amcl),
               0.7062363,
               tolerance = 0.001)
  expect_equal(septic_patients %>% portion_S(amcl, gent),
               0.9210074,
               tolerance = 0.001)
  expect_equal(septic_patients %>% portion_S(amcl, gent, also_single_tested = TRUE),
               0.9239669,
               tolerance = 0.001)

  # amcl+genta susceptibility around 92.1%
  expect_equal(suppressWarnings(rsi(septic_patients$amcl,
                                    septic_patients$gent,
                                    interpretation = "S")),
               0.9210074,
               tolerance = 0.000001)

  # percentages
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(R = portion_R(cipr, as_percent = TRUE),
                           I = portion_I(cipr, as_percent = TRUE),
                           S = portion_S(cipr, as_percent = TRUE),
                           n = n_rsi(cipr),
                           total = n()) %>%
                 pull(n) %>%
                 sum(),
               1409)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_p = portion_S(cipr, as_percent = TRUE),
                           cipro_n = n_rsi(cipr),
                           genta_p = portion_S(gent, as_percent = TRUE),
                           genta_n = n_rsi(gent),
                           combination_p = portion_S(cipr, gent, as_percent = TRUE),
                           combination_n = n_rsi(cipr, gent)) %>%
                 pull(combination_n),
               c(202, 488, 201, 499))

  expect_warning(portion_R(as.character(septic_patients$amcl)))
  expect_warning(portion_S(as.character(septic_patients$amcl)))
  expect_warning(portion_S(as.character(septic_patients$amcl,
                                             septic_patients$gent)))
  expect_warning(n_rsi(as.character(septic_patients$amcl,
                                    septic_patients$gent)))
  expect_equal(suppressWarnings(n_rsi(as.character(septic_patients$amcl,
                                                   septic_patients$gent))),
               1828)

  # check for errors
  expect_error(portion_IR("test", minimum = "test"))
  expect_error(portion_IR("test", as_percent = "test"))
  expect_error(portion_I("test", minimum = "test"))
  expect_error(portion_I("test", as_percent = "test"))
  expect_error(portion_S("test", minimum = "test"))
  expect_error(portion_S("test", as_percent = "test"))
  expect_error(portion_S("test", also_single_tested = "test"))

  # check too low amount of isolates
  expect_identical(suppressWarnings(portion_R(septic_patients$amox, minimum = nrow(septic_patients) + 1)),
                   NA)
  expect_identical(suppressWarnings(portion_I(septic_patients$amox, minimum = nrow(septic_patients) + 1)),
                   NA)
  expect_identical(suppressWarnings(portion_S(septic_patients$amox, minimum = nrow(septic_patients) + 1)),
                   NA)

  # warning for speed loss
  expect_warning(portion_R(as.character(septic_patients$gent)))
  expect_warning(portion_I(as.character(septic_patients$gent)))
  expect_warning(portion_S(septic_patients$amcl, as.character(septic_patients$gent)))

})

test_that("old rsi works", {
  # amox resistance in `septic_patients` should be around 58.53%
  expect_equal(suppressWarnings(rsi(septic_patients$amox)), 0.5853, tolerance = 0.0001)
  expect_equal(suppressWarnings(rsi(septic_patients$amox, interpretation = "S")), 1 - 0.5853, tolerance = 0.0001)

  # pita+genta susceptibility around 98.09%
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S",
                                    info = TRUE)),
               0.9498886,
               tolerance = 0.0001)

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
               c(202, 488, 201, 499))

  # portion_df
  expect_equal(
    septic_patients %>% select(amox) %>% portion_df() %>% pull(Value),
    c(septic_patients$amox %>% portion_S(),
      septic_patients$amox %>% portion_I(),
      septic_patients$amox %>% portion_R())
  )
  expect_equal(
    septic_patients %>% select(amox) %>% portion_df(combine_IR = TRUE) %>% pull(Value),
    c(septic_patients$amox %>% portion_S(),
      septic_patients$amox %>% portion_IR())
  )


})

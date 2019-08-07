# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("portion.R")

test_that("portions works", {
  # AMX resistance in `septic_patients`
  expect_equal(portion_R(septic_patients$AMX), 0.5557364, tolerance = 0.0001)
  expect_equal(portion_I(septic_patients$AMX), 0.002441009, tolerance = 0.0001)
  expect_equal(1 - portion_R(septic_patients$AMX) - portion_I(septic_patients$AMX),
               portion_S(septic_patients$AMX))
  expect_equal(portion_R(septic_patients$AMX) + portion_I(septic_patients$AMX),
               portion_IR(septic_patients$AMX))
  expect_equal(portion_S(septic_patients$AMX) + portion_I(septic_patients$AMX),
               portion_SI(septic_patients$AMX))

  expect_equal(septic_patients %>% portion_SI(AMC),
               0.7626397,
               tolerance = 0.0001)
  expect_equal(septic_patients %>% portion_SI(AMC, GEN),
               0.9408,
               tolerance = 0.0001)
  expect_equal(septic_patients %>% portion_SI(AMC, GEN, only_all_tested = TRUE),
               0.9382647,
               tolerance = 0.0001)

  # percentages
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(R = portion_R(CIP, as_percent = TRUE),
                           I = portion_I(CIP, as_percent = TRUE),
                           S = portion_S(CIP, as_percent = TRUE),
                           n = n_rsi(CIP),
                           total = n()) %>%
                 pull(n) %>%
                 sum(),
               1409)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_p = portion_SI(CIP, as_percent = TRUE),
                           cipro_n = n_rsi(CIP),
                           genta_p = portion_SI(GEN, as_percent = TRUE),
                           genta_n = n_rsi(GEN),
                           combination_p = portion_SI(CIP, GEN, as_percent = TRUE),
                           combination_n = n_rsi(CIP, GEN)) %>%
                 pull(combination_n),
               c(305, 617, 241, 711))

  expect_warning(portion_R(as.character(septic_patients$AMC)))
  expect_warning(portion_S(as.character(septic_patients$AMC)))
  expect_warning(portion_S(as.character(septic_patients$AMC,
                                        septic_patients$GEN)))
  expect_warning(n_rsi(as.character(septic_patients$AMC,
                                    septic_patients$GEN)))
  expect_equal(suppressWarnings(n_rsi(as.character(septic_patients$AMC,
                                                   septic_patients$GEN))),
               1879)

  # check for errors
  expect_error(portion_IR("test", minimum = "test"))
  expect_error(portion_IR("test", as_percent = "test"))
  expect_error(portion_I("test", minimum = "test"))
  expect_error(portion_I("test", as_percent = "test"))
  expect_error(portion_S("test", minimum = "test"))
  expect_error(portion_S("test", as_percent = "test"))
  expect_error(portion_S("test", also_single_tested = TRUE))

  # check too low amount of isolates
  expect_identical(suppressWarnings(portion_R(septic_patients$AMX, minimum = nrow(septic_patients) + 1)),
                   NA)
  expect_identical(suppressWarnings(portion_I(septic_patients$AMX, minimum = nrow(septic_patients) + 1)),
                   NA)
  expect_identical(suppressWarnings(portion_S(septic_patients$AMX, minimum = nrow(septic_patients) + 1)),
                   NA)

  # warning for speed loss
  expect_warning(portion_R(as.character(septic_patients$GEN)))
  expect_warning(portion_I(as.character(septic_patients$GEN)))
  expect_warning(portion_S(septic_patients$AMC, as.character(septic_patients$GEN)))

  # portion_df
  expect_equal(
    septic_patients %>% select(AMX) %>% portion_df() %>% pull(value),
    c(septic_patients$AMX %>% portion_SI(),
      septic_patients$AMX %>% portion_R())
  )
  expect_equal(
    septic_patients %>% select(AMX) %>% portion_df(combine_IR = TRUE) %>% pull(value),
    c(septic_patients$AMX %>% portion_S(),
      septic_patients$AMX %>% portion_IR())
  )
  expect_equal(
    septic_patients %>% select(AMX) %>% portion_df(combine_SI = FALSE) %>% pull(value),
    c(septic_patients$AMX %>% portion_S(),
      septic_patients$AMX %>% portion_I(),
      septic_patients$AMX %>% portion_R())
  )
  
  expect_error(portion_df(c("A", "B", "C")))
  expect_error(portion_df(septic_patients[,"date"]))
})

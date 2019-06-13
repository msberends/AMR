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

context("count.R")

test_that("counts work", {
  # AMX resistance in `septic_patients`
  expect_equal(count_R(septic_patients$AMX), 683)
  expect_equal(count_I(septic_patients$AMX), 3)
  expect_equal(count_S(septic_patients$AMX), 543)
  expect_equal(count_R(septic_patients$AMX) + count_I(septic_patients$AMX),
               count_IR(septic_patients$AMX))
  expect_equal(count_S(septic_patients$AMX) + count_I(septic_patients$AMX),
               count_SI(septic_patients$AMX))

  library(dplyr)
  expect_equal(septic_patients %>% count_S(AMC), 1342)
  expect_equal(septic_patients %>% count_S(AMC, GEN), 1660)
  expect_equal(septic_patients %>% count_all(AMC, GEN), 1798)
  expect_identical(septic_patients %>% count_all(AMC, GEN),
                   septic_patients %>% count_S(AMC, GEN) +
                     septic_patients %>% count_IR(AMC, GEN))

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro = count_S(CIP),
                           genta = count_S(GEN),
                           combination = count_S(CIP, GEN)) %>%
                 pull(combination),
               c(192, 446, 184, 474))

  # count_df
  expect_equal(
    septic_patients %>% select(AMX) %>% count_df() %>% pull(value),
    c(septic_patients$AMX %>% count_SI(),
      septic_patients$AMX %>% count_R())
  )
  expect_equal(
    septic_patients %>% select(AMX) %>% count_df(combine_IR = TRUE) %>% pull(value),
    c(septic_patients$AMX %>% count_S(),
      septic_patients$AMX %>% count_IR())
  )
  expect_equal(
    septic_patients %>% select(AMX) %>% count_df(combine_SI = FALSE) %>% pull(value),
    c(septic_patients$AMX %>% count_S(),
      septic_patients$AMX %>% count_I(),
      septic_patients$AMX %>% count_R())
  )

  # warning for speed loss
  expect_warning(count_R(as.character(septic_patients$AMC)))
  expect_warning(count_I(as.character(septic_patients$AMC)))
  expect_warning(count_S(as.character(septic_patients$AMC,
                                      septic_patients$GEN)))
  expect_warning(count_S(septic_patients$AMC,
                         as.character(septic_patients$GEN)))

  # check for errors
  expect_error(count_IR("test", minimum = "test"))
  expect_error(count_IR("test", as_percent = "test"))
  expect_error(count_I("test", minimum = "test"))
  expect_error(count_I("test", as_percent = "test"))
  expect_error(count_S("test", minimum = "test"))
  expect_error(count_S("test", as_percent = "test"))

  expect_error(count_df(c("A", "B", "C")))
  expect_error(count_df(septic_patients[,"date"]))

})

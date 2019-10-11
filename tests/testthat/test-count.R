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
  # AMX resistance in `example_isolates`
  expect_equal(count_R(example_isolates$AMX), 683)
  expect_equal(count_I(example_isolates$AMX), 3)
  expect_equal(count_S(example_isolates$AMX), 543)
  expect_equal(count_R(example_isolates$AMX) + count_I(example_isolates$AMX),
               count_IR(example_isolates$AMX))
  expect_equal(count_S(example_isolates$AMX) + count_I(example_isolates$AMX),
               count_SI(example_isolates$AMX))

  library(dplyr)
  expect_equal(example_isolates %>% count_S(AMC), 1342)
  expect_equal(example_isolates %>% count_S(AMC, GEN, only_all_tested = TRUE), 1660)
  expect_equal(example_isolates %>% count_S(AMC, GEN, only_all_tested = FALSE), 1728)
  expect_equal(example_isolates %>% count_all(AMC, GEN, only_all_tested = TRUE), 1798)
  expect_equal(example_isolates %>% count_all(AMC, GEN, only_all_tested = FALSE), 1936)
  expect_identical(example_isolates %>% count_all(AMC, GEN, only_all_tested = TRUE),
                   example_isolates %>% count_S(AMC, GEN, only_all_tested = TRUE) +
                     example_isolates %>% count_IR(AMC, GEN, only_all_tested = TRUE))

  # count of cases
  expect_equal(example_isolates %>%
                 group_by(hospital_id) %>%
                 summarise(cipro = count_SI(CIP),
                           genta = count_SI(GEN),
                           combination = count_SI(CIP, GEN)) %>%
                 pull(combination),
               c(253, 465, 192, 558))

  # count_df
  expect_equal(
    example_isolates %>% select(AMX) %>% count_df() %>% pull(value),
    c(example_isolates$AMX %>% count_SI(),
      example_isolates$AMX %>% count_R())
  )
  expect_equal(
    example_isolates %>% select(AMX) %>% count_df(combine_IR = TRUE) %>% pull(value),
    c(example_isolates$AMX %>% count_S(),
      example_isolates$AMX %>% count_IR())
  )
  expect_equal(
    example_isolates %>% select(AMX) %>% count_df(combine_SI = FALSE) %>% pull(value),
    c(example_isolates$AMX %>% count_S(),
      example_isolates$AMX %>% count_I(),
      example_isolates$AMX %>% count_R())
  )

  # warning for speed loss
  expect_warning(count_R(as.character(example_isolates$AMC)))
  expect_warning(count_I(as.character(example_isolates$AMC)))
  expect_warning(count_S(as.character(example_isolates$AMC,
                                      example_isolates$GEN)))
  expect_warning(count_S(example_isolates$AMC,
                         as.character(example_isolates$GEN)))

  # check for errors
  expect_error(count_IR("test", minimum = "test"))
  expect_error(count_IR("test", as_percent = "test"))
  expect_error(count_I("test", minimum = "test"))
  expect_error(count_I("test", as_percent = "test"))
  expect_error(count_S("test", minimum = "test"))
  expect_error(count_S("test", as_percent = "test"))

  expect_error(count_df(c("A", "B", "C")))
  expect_error(count_df(example_isolates[, "date"]))

})

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("first_isolate.R")

test_that("first isolates work", {
  
  skip_on_cran()
  
  # first isolates
  expect_equal(
    sum(
      first_isolate(x = example_isolates,
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    info = TRUE),
      na.rm = TRUE),
    1317)

  # first weighted isolates
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(x = example_isolates %>% mutate(keyab = key_antibiotics(.)),
                      # let syntax determine arguments automatically
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1413)
  # should be same for tibbles
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(x = example_isolates %>% dplyr::as_tibble() %>% mutate(keyab = key_antibiotics(.)),
                      # let syntax determine these automatically:
                      # col_date = "date",
                      # col_patient_id = "patient_id",
                      # col_mo = "mo",
                      # col_keyantibiotics = "keyab",
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1413)
  # when not ignoring I
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(x = example_isolates %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_mo = "mo",
                      col_keyantibiotics = "keyab",
                      ignore_I = FALSE,
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1436)
  # when using points
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(x = example_isolates %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_mo = "mo",
                      col_keyantibiotics = "keyab",
                      type = "points",
                      info = TRUE),
        na.rm = TRUE)),
    1417)

  # first non-ICU isolates
  expect_equal(
    sum(
      first_isolate(example_isolates,
                    col_mo = "mo",
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_icu = "ward_icu",
                    info = TRUE,
                    icu_exclude = TRUE),
      na.rm = TRUE),
    906)

  # set 1500 random observations to be of specimen type 'Urine'
  random_rows <- sample(x = 1:2000, size = 1500, replace = FALSE)
  expect_lt(
    sum(
      first_isolate(x = mutate(example_isolates,
                                 specimen = if_else(row_number() %in% random_rows,
                                                    "Urine",
                                                    "Other")),
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    col_specimen = "specimen",
                    filter_specimen = "Urine",
                    info = TRUE),
      na.rm = TRUE),
    1501)
  # same, but now exclude ICU
  expect_lt(
    sum(
      first_isolate(x = mutate(example_isolates,
                                 specimen = if_else(row_number() %in% random_rows,
                                                    "Urine",
                                                    "Other")),
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    col_specimen = "specimen",
                    filter_specimen = "Urine",
                    col_icu = "ward_icu",
                    icu_exclude = TRUE,
                    info = TRUE),
      na.rm = TRUE),
    1501)

  # "No isolates found"
  expect_message(example_isolates %>%
                   mutate(specimen = "test") %>%
                   mutate(first = first_isolate(., "date", "patient_id",
                                                col_mo = "mo",
                                                col_specimen = "specimen",
                                                filter_specimen = "something_unexisting",
                                                info = TRUE)))

  # printing of exclusion message
  expect_message(example_isolates %>%
                  first_isolate(col_date = "date",
                                col_mo = "mo",
                                col_patient_id = "patient_id",
                                col_testcode = "gender",
                                testcodes_exclude = "M",
                                info = TRUE))

  # errors
  expect_error(first_isolate("date", "patient_id", col_mo = "mo"))
  expect_error(first_isolate(example_isolates,
                             col_date = "non-existing col",
                             col_mo = "mo"))

  # look for columns itself
  expect_message(first_isolate(example_isolates))
  expect_message(first_isolate(example_isolates %>%
                               mutate(mo = as.character(mo)) %>%
                               left_join_microorganisms()))

  # if mo is not an mo class, result should be the same
  expect_identical(example_isolates %>%
                     mutate(mo = as.character(mo)) %>%
                     first_isolate(col_date = "date",
                                   col_mo = "mo",
                                   col_patient_id = "patient_id"),
                   example_isolates %>%
                     first_isolate(col_date = "date",
                                   col_mo = "mo",
                                   col_patient_id = "patient_id"))

  # missing dates should be no problem
  df <- example_isolates
  df[1:100, "date"] <- NA
  expect_equal(
    sum(
      first_isolate(x = df,
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    info = TRUE),
      na.rm = TRUE),
    1322)
  
  # unknown MOs
  expect_equal(example_isolates %>% 
                 mutate(mo = ifelse(mo == "B_ESCHR_COLI", "UNKNOWN", mo)) %>% 
                 mutate(first = first_isolate(., include_unknown = FALSE)) %>%
                 .$first %>% 
                 sum(),
               1062)
  expect_equal(example_isolates %>% 
                 mutate(mo = ifelse(mo == "B_ESCHR_COLI", "UNKNOWN", mo)) %>% 
                 mutate(first = first_isolate(., include_unknown = TRUE)) %>%
                 .$first %>% 
                 sum(),
               1529)
  expect_equal(example_isolates %>% 
                 mutate(mo = ifelse(mo == "B_ESCHR_COLI", NA, mo)) %>% 
                 mutate(first = first_isolate(.)) %>%
                 .$first %>% 
                 sum(),
               1062)

})

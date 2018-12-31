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

context("first_isolate.R")

test_that("first isolates work", {
  # septic_patients contains 1317 out of 2000 first isolates
  expect_equal(
    sum(
      first_isolate(tbl = septic_patients,
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    info = TRUE),
      na.rm = TRUE),
    1317)

  # septic_patients contains 1413 out of 2000 first *weighted* isolates
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      # let syntax determine these automatically:
                      # col_date = "date",
                      # col_patient_id = "patient_id",
                      # col_mo = "mo",
                      # col_keyantibiotics = "keyab",
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1413)
  # should be same for tibbles
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% dplyr::as_tibble() %>% mutate(keyab = key_antibiotics(.)),
                      # let syntax determine these automatically:
                      # col_date = "date",
                      # col_patient_id = "patient_id",
                      # col_mo = "mo",
                      # col_keyantibiotics = "keyab",
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1413)
  # and 1436 when not ignoring I
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_mo = "mo",
                      col_keyantibiotics = "keyab",
                      ignore_I = FALSE,
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1436)
  # and 1417 when using points
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_mo = "mo",
                      col_keyantibiotics = "keyab",
                      type = "points",
                      info = TRUE),
        na.rm = TRUE)),
    1417)

  # septic_patients contains 1163 out of 2000 first non-ICU isolates
  expect_equal(
    sum(
      first_isolate(septic_patients,
                    col_mo = "mo",
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_icu = "ward_icu",
                    info = TRUE,
                    icu_exclude = TRUE),
      na.rm = TRUE),
    1163)

  # set 1500 random observations to be of specimen type 'Urine'
  random_rows <- sample(x = 1:2000, size = 1500, replace = FALSE)
  expect_lt(
    sum(
      first_isolate(tbl = mutate(septic_patients,
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
      first_isolate(tbl = mutate(septic_patients,
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
  expect_message(septic_patients %>%
                   mutate(specimen = "test") %>%
                   mutate(first = first_isolate(., "date", "patient_id",
                                                col_mo = "mo",
                                                col_specimen = "specimen",
                                                filter_specimen = "something_unexisting")))

  # printing of exclusion message
  expect_output(septic_patients %>%
                  first_isolate(col_date = "date",
                                col_mo = "mo",
                                col_patient_id = "patient_id",
                                col_testcode = "gender",
                                testcodes_exclude = "M"))

  # errors
  expect_error(first_isolate("date", "patient_id", col_mo = "mo"))
  expect_error(first_isolate(septic_patients,
                             col_date = "non-existing col",
                             col_mo = "mo"))

  # look for columns itself
  expect_message(first_isolate(septic_patients))
  expect_error(first_isolate(septic_patients %>%
                               mutate(mo = as.character(mo)) %>%
                               left_join_microorganisms()))

  # if mo is not an mo class, result should be the same
  expect_identical(septic_patients %>%
                     mutate(mo = as.character(mo)) %>%
                     first_isolate(col_date = "date",
                                   col_mo = "mo",
                                   col_patient_id = "patient_id"),
                   septic_patients %>%
                     first_isolate(col_date = "date",
                                   col_mo = "mo",
                                   col_patient_id = "patient_id"))

})

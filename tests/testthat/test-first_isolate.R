context("first_isolate.R")

test_that("first isolates work", {
  # septic_patients contains 1331 out of 2000 first isolates
  expect_equal(
    sum(
      first_isolate(tbl = septic_patients,
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_mo = "mo",
                    info = TRUE),
      na.rm = TRUE),
    1331)

  # septic_patients contains 1426 out of 2000 first *weighted* isolates
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_mo = "mo",
                      col_keyantibiotics = "keyab",
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1426)
  # and 1449 when not ignoring I
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
    1449)
  # and 1430 when using points
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
    1430)

  # septic_patients contains 1176 out of 2000 first non-ICU isolates
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
    1176)

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
                                                filter_specimen = "something_unexisting",
                                                output_logical = FALSE)))

  # printing of exclusion message
  expect_output(septic_patients %>%
                            first_isolate(col_date = "date",
                                          col_mo = "mo",
                                          col_patient_id = "patient_id",
                                          col_testcode = "sex",
                                          testcodes_exclude = "M"))

  # errors
  expect_error(first_isolate("date", "patient_id", col_mo = "mo"))
  expect_error(first_isolate(septic_patients))
  expect_error(first_isolate(septic_patients,
                             col_date = "non-existing col",
                             col_mo = "mo"))

  expect_warning(septic_patients %>%
                   mutate(mo = as.character(mo)) %>%
                   first_isolate(col_date = "date",
                                 col_mo = "mo",
                                 col_patient_id = "patient_id"))

})

context("first_isolate.R")

test_that("first isolates work", {
  # septic_patients contains 1959 out of 2000 first isolates
  expect_equal(
    sum(
      first_isolate(tbl = septic_patients,
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_bactid = "bactid",
                    info = TRUE),
      na.rm = TRUE),
    1326)

  # septic_patients contains 1962 out of 2000 first *weighted* isolates
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_bactid = "bactid",
                      col_keyantibiotics = "keyab",
                      type = "keyantibiotics",
                      info = TRUE),
        na.rm = TRUE)),
    1421)
  # and 1961 when using points
  expect_equal(
    suppressWarnings(
      sum(
        first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                      col_date = "date",
                      col_patient_id = "patient_id",
                      col_bactid = "bactid",
                      col_keyantibiotics = "keyab",
                      type = "points",
                      info = TRUE),
        na.rm = TRUE)),
    1425)

  # septic_patients contains 1732 out of 2000 first non-ICU isolates
  expect_equal(
    sum(
      first_isolate(septic_patients,
                    col_bactid = "bactid",
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_icu = "ward_icu",
                    info = TRUE,
                    icu_exclude = TRUE),
      na.rm = TRUE),
    1171)

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
                    col_bactid = "bactid",
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
                    col_bactid = "bactid",
                    col_specimen = "specimen",
                    filter_specimen = "Urine",
                    col_icu = "ward_icu",
                    icu_exclude = TRUE,
                    info = TRUE),
      na.rm = TRUE),
    1501)
})

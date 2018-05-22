context("first_isolates.R")

test_that("keyantibiotics work", {
  expect_equal(length(key_antibiotics(septic_patients, info = FALSE)), nrow(septic_patients))
  expect_true(key_antibiotics_equal("SSS", "SSS"))
  expect_true(key_antibiotics_equal("SSS", "SIS", ignore_I = TRUE))
  expect_false(key_antibiotics_equal("SSS", "SIS", ignore_I = FALSE))
})

test_that("first isolates work", {
  # septic_patients contains 1960 out of 2000 first isolates
  #septic_ptns <- septic_patients
  expect_equal(sum(first_isolate(tbl = septic_patients,
                                 col_date = "date",
                                 col_patient_id = "patient_id",
                                 col_bactid = "bactid",
                                 info = FALSE)), 1960)

  # septic_patients contains 1962 out of 2000 first *weighted* isolates
  expect_equal(
    suppressWarnings(sum(
      first_isolate(tbl = septic_patients %>% mutate(keyab = key_antibiotics(.)),
                    col_date = "date",
                    col_patient_id = "patient_id",
                    col_bactid = "bactid",
                    col_keyantibiotics = "keyab",
                    type = "keyantibiotics",
                    info = TRUE))),
    1962)

  # septic_patients contains 1733 out of 2000 first non-ICU isolates
  expect_equal(
    sum(
      first_isolate(septic_patients, col_bactid = "bactid", col_date = "date", col_patient_id = "patient_id", col_icu = "ward_icu", info = TRUE, icu_exclude = TRUE)),
    1733
  )

  # set 1500 random observations to be of specimen type 'Urine'
  random_rows <- sample(x = 1:2000, size = 1500, replace = FALSE)
  expect_lt(sum(
    first_isolate(tbl = mutate(septic_patients,
                               specimen = if_else(row_number() %in% random_rows,
                                                  "Urine",
                                                  "Unknown")),
                  col_date = "date",
                  col_patient_id = "patient_id",
                  col_bactid = "bactid",
                  col_specimen = "specimen",
                  filter_specimen = "Urine",
                  info = TRUE)),
    1501)

})

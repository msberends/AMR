context("key_antibiotics.R")

test_that("keyantibiotics work", {
  expect_equal(length(key_antibiotics(septic_patients, warnings = FALSE)), nrow(septic_patients))
  expect_false(all(is.na(key_antibiotics(septic_patients))))
  expect_true(key_antibiotics_equal("SSS", "SSS"))
  expect_false(key_antibiotics_equal("SSS", "SRS"))
  expect_true(key_antibiotics_equal("SSS", "SIS", ignore_I = TRUE))
  expect_false(key_antibiotics_equal("SSS", "SIS", ignore_I = FALSE))
  expect_true(key_antibiotics_equal(".SS", "SI.", ignore_I = TRUE))
  expect_false(key_antibiotics_equal(".SS", "SI.", ignore_I = FALSE))
})

context("ab_property.R")

test_that("ab_property works", {
  expect_equal(ab_certe("amox"), "amox")
  expect_equal(ab_official("amox"), "Amoxicillin")
  expect_equal(ab_official_nl("amox"), "Amoxicilline")
  expect_equal(ab_trivial_nl("amox"), "Amoxicilline")
  expect_equal(ab_umcg("amox"), "AMOX")
})

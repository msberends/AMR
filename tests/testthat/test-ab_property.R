context("ab_property.R")

test_that("ab_property works", {
  expect_equal(ab_certe("amox"), "amox")
  expect_equal(ab_name("amox", language = "en"), "Amoxicillin")
  expect_equal(ab_name("amox", language = "nl"), "Amoxicilline")
  expect_equal(ab_official("amox", language = "en"), "Amoxicillin")
  expect_equal(ab_trivial_nl("amox"), "Amoxicilline")
  expect_equal(ab_umcg("amox"), "AMOX")
  expect_equal(class(ab_tradenames("amox")), "character")
  expect_equal(class(ab_tradenames(c("amox", "amox"))), "list")
  expect_equal(ab_atc("amox"), as.character(as.atc("amox")))

  expect_error(ab_property("amox", "invalid property"))
  expect_error(ab_name("amox", language = "INVALID"))
  expect_output(print(ab_name("amox", language = NULL)))
})

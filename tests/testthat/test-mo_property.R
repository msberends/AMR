context("mo_property.R")

test_that("mo_property works", {
  expect_equal(mo_family("E. coli"), "Enterobacteriaceae")
  expect_equal(mo_genus("E. coli"), "Escherichia")
  expect_equal(mo_species("E. coli"), "coli")
  expect_equal(mo_subspecies("E. coli"), NA_character_)
  expect_equal(mo_fullname("E. coli"), "Escherichia coli")
  expect_equal(mo_type("E. coli"), "Bacteria")
  expect_equal(mo_gramstain("E. coli"), "Negative rods")
  expect_equal(mo_aerobic("E. coli"), TRUE)
  expect_equal(mo_type_nl("E. coli"), "Bacterie")
  expect_equal(mo_gramstain_nl("E. coli"), "Negatieve staven")
})

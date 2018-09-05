context("mo_property.R")

test_that("mo_property works", {
  expect_equal(mo_family("E. coli"), "Enterobacteriaceae")
  expect_equal(mo_genus("E. coli"), "Escherichia")
  expect_equal(mo_species("E. coli"), "coli")
  expect_equal(mo_subspecies("E. coli"), "")
  expect_equal(mo_fullname("E. coli"), "Escherichia coli")
  expect_equal(mo_type("E. coli"), "Bacteria")
  expect_equal(mo_gramstain("E. coli"), "Negative rods")
  expect_equal(mo_aerobic("E. coli"), TRUE)

  expect_equal(mo_shortname("MRSA"), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = TRUE), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = "all"), "CoPS")
  expect_equal(mo_shortname("S. aga"), "S. agalactiae")
  expect_equal(mo_shortname("S. aga", Lancefield = TRUE), "GBS")

  expect_equal(mo_type("E. coli", language = "de"), "Bakterien")
  expect_equal(mo_gramstain("E. coli", language = "de"), "Negative Staebchen")

  expect_equal(mo_type("E. coli", language = "nl"), "Bacterie")
  expect_equal(mo_gramstain("E. coli", language = "nl"), "Negatieve staven")

  expect_error(mo_type("E. coli", language = "INVALID"))
})

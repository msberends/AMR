context("mo_property.R")

test_that("mo_property works", {
  expect_equal(mo_subkingdom("E. coli"), "Negibacteria")
  expect_equal(mo_phylum("E. coli"), "Proteobacteria")
  expect_equal(mo_class("E. coli"), "Gammaproteobacteria")
  expect_equal(mo_order("E. coli"), "Enterobacteriales")
  expect_equal(mo_family("E. coli"), "Enterobacteriaceae")
  expect_equal(mo_genus("E. coli"), "Escherichia")
  expect_equal(mo_species("E. coli"), "coli")
  expect_equal(mo_subspecies("E. coli"), "")
  expect_equal(mo_fullname("E. coli"), "Escherichia coli")
  expect_equal(mo_type("E. coli", language = "en"), "Bacteria")
  expect_equal(mo_gramstain("E. coli", language = "en"), "Gram negative")
  expect_equal(class(mo_taxonomy("E. coli")), "list")

  expect_equal(mo_shortname("MRSA"), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = TRUE), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = "all"), "CoPS")
  expect_equal(mo_shortname("S. aga"), "S. agalactiae")
  expect_equal(mo_shortname("S. aga", Lancefield = TRUE), "GBS")

  # test integrity
  # library(dplyr)
  # rnd <- sample(1:nrow(AMR::microorganisms), 500, replace = FALSE) # random 500 rows
  # MOs <- AMR::microorganisms %>% filter(!is.na(mo),
  #                                       species != "species",
  #                                       dplyr::row_number() %in% rnd)
  # expect_identical(MOs$fullname, mo_fullname(MOs$fullname, language = "en"))

  # check languages
  expect_equal(mo_type("E. coli", language = "de"), "Bakterien")
  expect_equal(mo_gramstain("E. coli", language = "nl"), "Gram-negatief")

  expect_output(print(mo_gramstain("E. coli", language = "en")))
  expect_output(print(mo_gramstain("E. coli", language = "de")))
  expect_output(print(mo_gramstain("E. coli", language = "nl")))
  expect_output(print(mo_gramstain("E. coli", language = "es")))
  expect_output(print(mo_gramstain("E. coli", language = "pt")))
  expect_output(print(mo_gramstain("E. coli", language = "it")))
  expect_output(print(mo_gramstain("E. coli", language = "fr")))

  expect_error(mo_gramstain("E. coli", language = "UNKNOWN"))

})

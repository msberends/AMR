context("mo_property.R")

test_that("mo_property works", {
  expect_equal(mo_family("E. coli"), "Enterobacteriaceae")
  expect_equal(mo_genus("E. coli"), "Escherichia")
  expect_equal(mo_species("E. coli"), "coli")
  expect_equal(mo_subspecies("E. coli"), "")
  expect_equal(mo_fullname("E. coli"), "Escherichia coli")
  expect_equal(mo_type("E. coli", language = "en"), "Bacteria")
  expect_equal(mo_gramstain("E. coli", language = "en"), "Negative rods")
  expect_equal(mo_aerobic("E. coli"), TRUE)

  expect_equal(mo_shortname("MRSA"), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = TRUE), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = "all"), "CoPS")
  expect_equal(mo_shortname("S. aga"), "S. agalactiae")
  expect_equal(mo_shortname("S. aga", Lancefield = TRUE), "GBS")

  # test integrity
  library(dplyr)
  MOs <- AMR::microorganisms %>% filter(!is.na(mo))
  expect_identical(MOs$fullname, mo_fullname(MOs$fullname, language = "en"))

  mo_clean <- MOs$mo
  mo_from_shortname <- as.mo(mo_shortname(mo_clean))
  mo_clean <- mo_clean[nchar(mo_from_shortname) == 6 &
                         !is.na(mo_from_shortname) &
                         !mo_from_shortname %like% "...SPP"]
  mo_from_shortname <- mo_from_shortname[nchar(mo_from_shortname) == 6 &
                                           !is.na(mo_from_shortname) &
                                           !mo_from_shortname %like% "...SPP"]
  tb <- tibble(a = substr(mo_clean, 1, 6),
               b = mo_from_shortname,
               c = a == b,
               d = mo_shortname(a),
               e = mo_shortname(b),
               f = d == e)
  expect_gt(sum(tb$c) / nrow(tb), 0.9) # more than 90% of MO code should be identical
  expect_identical(sum(tb$f), nrow(tb)) # all shortnames should be identical

  # check languages
  expect_equal(mo_type("E. coli", language = "de"), "Bakterium")
  expect_equal(mo_type("E. coli", language = "nl"), "Bacterie")
  expect_equal(mo_gramstain("E. coli", language = "nl"), "Negatieve staven")

  expect_output(print(mo_gramstain("E. coli", language = "en")))
  expect_output(print(mo_gramstain("E. coli", language = "de")))
  expect_output(print(mo_gramstain("E. coli", language = "nl")))
  expect_output(print(mo_gramstain("E. coli", language = "es")))
  expect_output(print(mo_gramstain("E. coli", language = "pt")))
  expect_output(print(mo_gramstain("E. coli", language = "it")))
  expect_output(print(mo_gramstain("E. coli", language = "fr")))

  expect_error(mo_gramstain("E. coli", language = "UNKNOWN"))

})

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

context("mo_property.R")

test_that("mo_property works", {
  expect_equal(mo_kingdom("E. coli"), "Bacteria")
  expect_equal(mo_subkingdom("E. coli"), "Negibacteria")
  expect_equal(mo_phylum("E. coli"), "Proteobacteria")
  expect_equal(mo_class("E. coli"), "Gammaproteobacteria")
  expect_equal(mo_order("E. coli"), "Enterobacteriales")
  expect_equal(mo_family("E. coli"), "Enterobacteriaceae")
  expect_equal(mo_genus("E. coli"), "Escherichia")
  expect_equal(mo_species("E. coli"), "coli")
  expect_equal(mo_subspecies("E. coli"), NA_character_)
  expect_equal(mo_fullname("E. coli"), "Escherichia coli")
  expect_equal(mo_type("E. coli", language = "en"), "Bacteria")
  expect_equal(mo_gramstain("E. coli", language = "en"), "Gram negative")
  expect_equal(mo_TSN("E. coli"), 285)
  expect_equal(class(mo_taxonomy("E. coli")), "list")
  expect_equal(names(mo_taxonomy("E. coli")), c("kingdom", "subkingdom", "phylum", "class", "order",
                                                "family", "genus", "species", "subspecies"))

  expect_equal(mo_ref("E. coli"), "Castellani and Chalmers, 1919")
  expect_equal(mo_authors("E. coli"), "Castellani and Chalmers")
  expect_equal(mo_year("E. coli"), 1919)

  expect_equal(mo_shortname("MRSA"), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = TRUE), "S. aureus")
  expect_equal(mo_shortname("MRSA", Becker = "all", language = "en"), "CoPS")
  expect_equal(mo_shortname("S. aga"), "S. agalactiae")
  expect_equal(mo_shortname("S. aga", Lancefield = TRUE), "GBS")

  # test integrity
  MOs <- AMR::microorganisms
  expect_identical(MOs$fullname, mo_fullname(MOs$fullname, language = "en"))

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

  # manual property function
  expect_error(mo_property("E. coli", property = c("tsn", "fullname")))
  expect_error(mo_property("E. coli", property = "UNKNOWN"))
  expect_identical(mo_property("E. coli", property = "fullname"),
                   mo_fullname("E. coli"))
  expect_identical(mo_property("E. coli", property = "genus"),
                   mo_genus("E. coli"))
  expect_identical(mo_property("E. coli", property = "species"),
                   mo_species("E. coli"))

  # check vector with random values
  library(dplyr)
  df_sample <- AMR::microorganisms %>% sample_n(100)
  expect_identical(df_sample %>% pull(mo) %>% mo_fullname(language = "en"),
                   df_sample %>% pull(fullname))

})

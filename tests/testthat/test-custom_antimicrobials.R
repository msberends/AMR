# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in the Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

test_that("test-custom ab.R", {
  skip_on_cran()

  ab_reset_session()

  expect_message(as.ab("testab", info = TRUE))

  suppressMessages(
    add_custom_antimicrobials(
      data.frame(
        ab = "TESTAB",
        name = "Test Antibiotic",
        group = "Test Group"
      )
    )
  )

  expect_identical(as.character(as.ab("testab")), "TESTAB")
  expect_identical(ab_name("testab"), "Test Antibiotic")
  expect_identical(ab_group("testab"), "Test Group")
})

test_that("test-custom ab synonyms", {
  skip_on_cran()

  suppressMessages(clear_custom_antimicrobials())
  suppressMessages(ab_reset_session())

  # codes must exist
  expect_error(add_custom_antimicrobial_synonyms("NOTANAB", "something"))

  # ASCII and non-ASCII synonyms (a Korean trade name of piperacillin/tazobactam)
  tazocin_ko <- "\uD0C0\uC870\uC2E0\uC8FC"
  suppressMessages(add_custom_antimicrobial_synonyms("TZP", c("Tazocinum Testname", tazocin_ko)))
  expect_identical(as.character(as.ab("Tazocinum Testname")), "TZP")
  expect_identical(as.character(as.ab(tazocin_ko)), "TZP")
  expect_identical(as.character(as.ab("\uD0C0\uC870 \uC2E0\uC8FC")), "TZP") # white space is ignored
  expect_identical(as.character(as.ab(c(tazocin_ko, "amoxicillin", NA))), c("TZP", "AMX", NA))
  expect_identical(ab_name(tazocin_ko), ab_name("TZP"))

  # another non-ASCII name that was not added must not be matched to TZP
  expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab("\uBA54\uB85C\uD39C")))), "TZP"))

  # data.frame input, also for custom antimicrobials
  suppressMessages(add_custom_antimicrobials(data.frame(ab = "TESTSYN", name = "Test Synonym Antibiotic")))
  suppressMessages(add_custom_antimicrobial_synonyms(data.frame(ab = "TESTSYN", synonym = "Testosyn")))
  expect_identical(as.character(as.ab("testosyn")), "TESTSYN")

  # one synonym cannot refer to two antimicrobials
  expect_error(add_custom_antimicrobial_synonyms("MEM", tazocin_ko))

  # clearing removes the synonyms as well
  suppressMessages(clear_custom_antimicrobials())
  expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab(tazocin_ko)))), "TZP"))
})

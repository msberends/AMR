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

test_that("test-italicise_taxonomy.R", {
  skip_on_cran()

  expect_identical(
    italicise_taxonomy("test for E. coli"),
    "test for *E. coli*"
  )
  expect_identical(
    italicise_taxonomy("test for E. coli"),
    italicize_taxonomy("test for E. coli")
  )
  if (AMR:::has_colour()) {
    expect_identical(
      italicise_taxonomy("test for E. coli", type = "ansi"),
      "test for \033[3mE. coli\033[23m"
    )
  }
  expect_identical(
    italicise_taxonomy("test for E. coli", "html"),
    "test for <i>E. coli</i>"
  )

  # species only after a genus or abbreviated genus, punctuation outside the italics
  expect_identical(
    italicise_taxonomy(c(
      "Escherichia coli", "coli bacteria", "(S. aureus)", "E. coli,", "Other",
      "Salmonella Typhimurium", "Streptococcus group A", "Staphylococcus spp.",
      "Klebsiella pneumoniae subsp. ozaenae", "K. pneumoniae ozaenae", "Enterobacteriaceae", "", NA
    )),
    c(
      "*Escherichia coli*", "coli bacteria", "(*S. aureus*)", "*E. coli*,", "Other",
      "*Salmonella* Typhimurium", "*Streptococcus* group A", "*Staphylococcus* spp.",
      "*Klebsiella pneumoniae* subsp. *ozaenae*", "*K. pneumoniae ozaenae*", "*Enterobacteriaceae*", "", NA
    )
  )

  # plotmath, also used by labels_taxonomy()
  expect_identical(
    italicise_taxonomy(c("E. coli", "Other", "test for E. coli", NA), type = "plotmath"),
    expression(italic("E. coli"), "Other", paste("test for ", italic("E. coli")), NA_character_)
  )
  expect_identical(
    labels_taxonomy(factor(c("E. coli", "Other"))),
    expression(italic("E. coli"), "Other")
  )
  expect_identical(labels_taxonomy(character(0)), expression())

  # used internally by the plotting functions
  expect_true(AMR:::contains_taxonomy(c("E. coli", "Other", NA)))
  expect_false(AMR:::contains_taxonomy(c("Other", "CoNS", NA)))
  expect_identical(
    AMR:::labeller_taxonomy(list(mo = c("E. coli", "CoNS"))),
    list(mo = list(expression(italic("E. coli")), expression("CoNS")))
  )
})

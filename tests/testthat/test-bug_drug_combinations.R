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

test_that("test-bug/drug.R", {
  skip_on_cran()

  b <- suppressWarnings(bug_drug_combinations(example_isolates))
  expect_inherits(b, "bug_drug_combinations")
  expect_output(suppressMessages(print(b)))
  expect_true(is.data.frame(format(b)))
  expect_true(is.data.frame(format(b, add_ab_group = FALSE)))

  # format() calculates resistance following the guideline, taken from bug_drug_combinations() unless set
  b_eucast <- suppressMessages(bug_drug_combinations(example_isolates, FUN = function(x) x, guideline = "EUCAST"))
  b_clsi <- suppressMessages(bug_drug_combinations(example_isolates, FUN = function(x) x, guideline = "CLSI"))
  expect_identical(attr(b_clsi, "guideline"), "CLSI")
  cons_gen <- subset(b_eucast, mo == "B_STPHY_CONS" & ab == "GEN")
  n_R <- cons_gen$R + cons_gen$NWT + cons_gen$NS
  n_IR <- n_R + cons_gen$I
  get_cons_gen <- function(f) f[f$Drug %like% "GEN", "B_STPHY_CONS", drop = TRUE]
  expect_true(grepl(paste0("(", n_R, "/", cons_gen$total, ")"), get_cons_gen(format(b_eucast, translate_ab = "ab", language = "en")), fixed = TRUE))
  expect_true(grepl(paste0("(", n_IR, "/", cons_gen$total, ")"), get_cons_gen(format(b_clsi, translate_ab = "ab", language = "en")), fixed = TRUE))
  expect_identical(format(b_eucast, guideline = "CLSI"), format(b_clsi))
  # with guideline = NULL, format() requires its own guideline
  b_none <- suppressMessages(bug_drug_combinations(example_isolates, FUN = function(x) x, guideline = NULL))
  expect_error(format(b_none))
  expect_identical(format(b_none, guideline = "CLSI"), format(b_clsi))
  # the deprecated combine_SI must be followed literally if set, regardless of the guideline
  expect_identical(suppressWarnings(format(b_clsi, combine_SI = TRUE)), format(b_eucast))
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_true(example_isolates %>%
      group_by(ward) %>%
      bug_drug_combinations(FUN = mo_gramstain) %>%
      is.data.frame())
  }
})

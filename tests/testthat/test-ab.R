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
# Center Groningen in The Netherlands, in collaboration with many      #
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

test_that("test-ab.R", {
  skip_on_cran()

  ab_reset_session()

  expect_equal(
    as.character(as.ab(c(
      "J01FA01",
      "J 01 FA 01",
      "Erythromycin",
      "eryt",
      "ERYT",
      "ERY",
      "erytromicine",
      "Erythrocin"
    ))),
    rep("ERY", 8)
  )

  expect_identical(class(as.ab("amox")), c("ab", "character"))
  expect_identical(class(AMR::antimicrobials$ab), c("ab", "character"))
  expect_true(is.ab(as.ab("amox")))
  expect_output(print(as.ab("amox")))
  expect_output(print(data.frame(a = as.ab("amox"))))

  expect_warning(as.ab("J00AA00")) # ATC not yet available in data set
  # expect_warning(as.ab("UNKNOWN"))

  expect_output(print(as.ab("amox")))

  expect_equal(
    as.character(as.ab("Phloxapen")),
    "FLC"
  )

  expect_equal(
    suppressWarnings(as.character(as.ab(c("Bacteria", "Bacterial")))),
    c(NA, "TMP")
  )

  expect_equal(
    as.character(as.ab("Amoxy + clavulaanzuur")),
    "AMC"
  )

  expect_equal(
    as.character(as.ab(c("mreopenem", "co-maoxiclav"))),
    c("MEM", "AMC")
  )

  # expect_warning(as.ab("cipro mero"))

  # based on Levenshtein distance
  expect_identical(ab_name("ceftazidim/avibactam", language = NULL), "Ceftazidime/avibactam")

  expect_identical(
    as.character(as.ab(c(
      "gentamicine High Level",
      "gentamicine High",
      "gentamicine (High Level)",
      "gentamicine (High)",
      "gentamicine HL",
      "gentamicine H-L",
      "gentamicine (HL)",
      "gentamicine (H-L)"
    ))),
    rep("GEH", 8)
  )

  # assigning and subsetting
  x <- AMR::antimicrobials$ab
  expect_inherits(x[1], "ab")
  expect_inherits(x[[1]], "ab")
  expect_inherits(c(x[1], x[9]), "ab")
  expect_inherits(unique(x[1], x[9]), "ab")
  expect_inherits(rep(x[1], 2), "ab")
  # expect_warning(x[1] <- "invalid code")
  # expect_warning(x[[1]] <- "invalid code")
  # expect_warning(c(x[1], "test"))
})

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

test_that("test-ab_property.R", {
  skip_on_cran()

  ab_reset_session()

  expect_identical(ab_name("AMX", language = NULL), "Amoxicillin")
  expect_identical(ab_atc("AMX"), c("J01CA04", "QG51AA03", "QJ01CA04"))
  expect_identical(ab_cid("AMX"), as.integer(33613))

  expect_inherits(ab_tradenames("AMX"), "character")
  expect_inherits(ab_tradenames(c("AMX", "AMX")), "list")

  expect_identical(ab_group("AMX", language = NULL), "Beta-lactams/penicillins")
  expect_identical(ab_atc_group1("AMX", language = NULL), "Beta-lactam antibacterials, penicillins")
  expect_identical(ab_atc_group2("AMX", language = NULL), "Penicillins with extended spectrum")

  expect_identical(ab_name("Fluclox", language = NULL), "Flucloxacillin")
  expect_identical(ab_name("fluklox", language = NULL), "Flucloxacillin")
  expect_identical(ab_name("floxapen", language = NULL), "Flucloxacillin")
  expect_identical(ab_name(21319, language = NULL), "Flucloxacillin")
  expect_identical(ab_name("J01CF05", language = NULL), "Flucloxacillin")

  expect_identical(ab_ddd("AMX", "oral"), 1.5)
  expect_identical(ab_ddd_units("AMX", "iv"), "g")
  expect_identical(ab_ddd("AMX", "iv"), 3)

  expect_identical(ab_name(x = c("AMC", "PLB"), language = NULL), c("Amoxicillin/clavulanic acid", "Polymyxin B"))
  expect_identical(
    ab_name(x = c("AMC", "PLB"), tolower = TRUE, language = NULL),
    c("amoxicillin/clavulanic acid", "polymyxin B")
  )

  expect_inherits(ab_info("AMX"), "list")

  expect_error(ab_property("amox", "invalid property"))
  expect_error(ab_name("amox", language = "INVALID"))
  expect_output(print(ab_name("amox", language = NULL)))

  expect_equal(ab_name("21066-6", language = NULL), "Ampicillin")
  expect_equal(
    ab_loinc("ampicillin"),
    c("101477-8", "101478-6", "18864-9", "18865-6", "20374-5", "21066-6", "23618-2", "27-3", "28-1", "29-9", "30-7", "31-5", "32-3", "33-1", "3355-5", "33562-0", "33919-2", "34-9", "43883-8", "43884-6", "6979-9", "6980-7", "87604-5")
  )

  expect_true(ab_url("AMX") %like% "fhi[.]no")

  expect_identical(
    colnames(set_ab_names(example_isolates[, 17:22])),
    c("cefoxitin", "cefotaxime", "ceftazidime", "ceftriaxone", "gentamicin", "tobramycin")
  )
  expect_identical(
    colnames(set_ab_names(example_isolates[, 17:22], language = "nl", snake_case = FALSE)),
    c("Cefoxitine", "Cefotaxim", "Ceftazidim", "Ceftriaxon", "Gentamicine", "Tobramycine")
  )
  expect_identical(
    colnames(set_ab_names(example_isolates[, 17:22], property = "atc")),
    c("J01DC01", "J01DD01", "J01DD02", "J01DD04", "J01GB03", "J01GB01")
  )

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_identical(
      example_isolates %>% set_ab_names(),
      example_isolates %>% rename_with(set_ab_names)
    )
    expect_true(all(c(
      "SXT", "nitrofurantoin", "fosfomycin", "linezolid", "ciprofloxacin",
      "moxifloxacin", "vancomycin", "TEC"
    ) %in%
      (example_isolates %>%
        set_ab_names(NIT:VAN) %>%
        colnames())))
  }
})

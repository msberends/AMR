# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

expect_equal(AMR:::percentage(0.25), "25%")
expect_equal(AMR:::percentage(0.5), "50%")
expect_equal(AMR:::percentage(0.500, digits = 1), "50.0%")
expect_equal(AMR:::percentage(0.1234), "12.3%")
# round up 0.5
expect_equal(AMR:::percentage(0.0054), "0.5%")
expect_equal(AMR:::percentage(0.0055), "0.6%")

# test functions on all R versions -  R < 3.3 did not contain these
expect_equal(strrep("A", 5), "AAAAA")
expect_equal(strrep(c("A", "B"), c(5, 2)), c("AAAAA", "BB"))
expect_equal(trimws(" test "), "test")
expect_equal(trimws(" test ", "l"), "test ")
expect_equal(trimws(" test ", "r"), " test")
expect_equal(trimws2(" test "), "test")
expect_equal(trimws2(" test ", "l"), "test ")
expect_equal(trimws2(" test ", "r"), " test")

expect_warning(AMR:::generate_warning_abs_missing(c("AMP", "AMX")))
expect_warning(AMR:::generate_warning_abs_missing(c("AMP", "AMX"), any = TRUE))
expect_warning(AMR:::get_column_abx(example_isolates, hard_dependencies = "FUS"))
expect_message(AMR:::get_column_abx(example_isolates, soft_dependencies = "FUS"))

if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  expect_warning(AMR:::get_column_abx(rename(example_isolates, thisone = AMX), amox = "thisone", tmp = "thisone", verbose = TRUE))
  expect_warning(AMR:::get_column_abx(rename(example_isolates, thisone = AMX), amox = "thisone", tmp = "thisone", verbose = FALSE))
}

# we rely on "grouped_tbl" being a class of grouped tibbles, so run a test that checks for this:
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  expect_true(AMR:::is_null_or_grouped_tbl(example_isolates %>% group_by(ward)))
}

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

test_that("test-misc.R", {
  skip_on_cran()

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
  expect_equal(AMR:::trimws2(" test "), "test")
  expect_equal(AMR:::trimws2(" test ", "l"), "test ")
  expect_equal(AMR:::trimws2(" test ", "r"), " test")

  expect_message(AMR:::get_column_abx(example_isolates, soft_dependencies = "FUS"))

  # we rely on "grouped_tbl" being a class of grouped tibbles, so run a test that checks for this:
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_true(AMR:::is_null_or_grouped_tbl(example_isolates %>% group_by(ward)))
  }


  # test get_current_data() ----

  is_right <- FALSE
  check_df <- function(check_element, return_val = 0) {
    is_right <<- FALSE
    for (env in sys.frames()) {
      if (!is.null(env[[check_element]]) && is.data.frame(env[[check_element]])) {
        is_right <<- all(colnames(example_isolates) %in% colnames(env[[check_element]]))
      }
    }
    return_val
  }

  df <- example_isolates[, check_df("x")]
  expect_true(is_right, info = "the environmental data cannot be found for base `x`")

  # should work on R >=3.6.3 or so
  df <- example_isolates[c(1:3), check_df("x")]
  if (!is_right) {
    # otherwise, this is needed for older versions
    df <- example_isolates[c(1:3), check_df("xx")]
    expect_true(is_right, info = "the environmental data cannot be found for base `x` or `xx`") 
  }

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    df <- example_isolates %>% select(mo, check_df("data123"))
    expect_false(is_right, info = "just a check if non-sense is not being gathered by get_current_data()")

    df <- example_isolates %>% select(mo, check_df(".data"))
    expect_true(is_right, info = "the environmental data cannot be found for dplyr/select()")

    df <- example_isolates %>% select_at(check_df(".tbl"))
    expect_true(is_right, info = "the environmental data cannot be found for dplyr/select_at()")
  }

  if (AMR:::pkg_is_available("tidymodels", also_load = TRUE)) {
    resistance_recipe <- recipe(mo ~ ., data = example_isolates) %>%
      step_corr(check_df("training")) %>%
      prep()
    expect_true(is_right, info = "the environmental data cannot be found for tidymodels/prep()")
  }
})

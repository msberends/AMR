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

test_that("test-proportion.R", {
  skip_on_cran()

  expect_equal(proportion_R(example_isolates$AMX), resistance(example_isolates$AMX))
  expect_equal(proportion_SI(example_isolates$AMX), susceptibility(example_isolates$AMX))
  # AMX resistance in `example_isolates`
  expect_equal(proportion_R(example_isolates$AMX), 0.5955556, tolerance = 0.0001)
  expect_equal(proportion_I(example_isolates$AMX), 0.002222222, tolerance = 0.0001)
  expect_equal(sir_confidence_interval(example_isolates$AMX)[1], 0.5688204, tolerance = 0.0001)
  expect_equal(sir_confidence_interval(example_isolates$AMX)[2], 0.6218738, tolerance = 0.0001)
  expect_equal(
    1 - proportion_R(example_isolates$AMX) - proportion_I(example_isolates$AMX),
    proportion_S(example_isolates$AMX)
  )
  expect_equal(
    proportion_R(example_isolates$AMX) + proportion_I(example_isolates$AMX),
    proportion_IR(example_isolates$AMX)
  )
  expect_equal(
    proportion_S(example_isolates$AMX) + proportion_I(example_isolates$AMX),
    proportion_SI(example_isolates$AMX)
  )

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_equal(example_isolates %>% proportion_SI(AMC),
      0.7626397,
      tolerance = 0.0001
    )
    expect_equal(example_isolates %>% proportion_SI(AMC, GEN),
      0.9408,
      tolerance = 0.0001
    )
    expect_equal(example_isolates %>% proportion_SI(AMC, GEN, only_all_tested = TRUE),
      0.9382647,
      tolerance = 0.0001
    )

    # percentages
    expect_equal(
      example_isolates %>%
        group_by(ward) %>%
        summarise(
          R = proportion_R(CIP, as_percent = TRUE),
          I = proportion_I(CIP, as_percent = TRUE),
          S = proportion_S(CIP, as_percent = TRUE),
          n = n_sir(CIP),
          total = n()
        ) %>%
        pull(n) %>%
        sum(),
      1409
    )

    # count of cases
    expect_equal(
      example_isolates %>%
        group_by(ward) %>%
        summarise(
          cipro_p = proportion_SI(CIP, as_percent = TRUE),
          cipro_n = n_sir(CIP),
          genta_p = proportion_SI(GEN, as_percent = TRUE),
          genta_n = n_sir(GEN),
          combination_p = proportion_SI(CIP, GEN, as_percent = TRUE),
          combination_n = n_sir(CIP, GEN)
        ) %>%
        pull(combination_n),
      c(1181, 577, 116)
    )

    # proportion_df
    expect_equal(
      example_isolates %>% select(AMX) %>% proportion_df() %>% pull(value),
      c(
        example_isolates$AMX %>% proportion_SI(),
        example_isolates$AMX %>% proportion_R()
      )
    )
    expect_equal(
      example_isolates %>% select(AMX) %>% proportion_df(combine_SI = FALSE) %>% pull(value),
      c(
        example_isolates$AMX %>% proportion_S(),
        example_isolates$AMX %>% proportion_I(),
        example_isolates$AMX %>% proportion_R()
      )
    )

    # expect_warning(example_isolates %>% group_by(ward) %>% summarise(across(KAN, sir_confidence_interval)))
  }

  # expect_warning(proportion_R(as.character(example_isolates$AMC)))
  # expect_warning(proportion_S(as.character(example_isolates$AMC)))
  # expect_warning(proportion_S(as.character(example_isolates$AMC, example_isolates$GEN)))

  # expect_warning(n_sir(as.character(example_isolates$AMC, example_isolates$GEN)))
  expect_equal(
    suppressWarnings(n_sir(as.character(
      example_isolates$AMC,
      example_isolates$GEN
    ))),
    1879
  )

  # check for errors
  expect_error(proportion_IR("test", minimum = "test"))
  expect_error(proportion_IR("test", as_percent = "test"))
  expect_error(proportion_I("test", minimum = "test"))
  expect_error(proportion_I("test", as_percent = "test"))
  expect_error(proportion_S("test", minimum = "test"))
  expect_error(proportion_S("test", as_percent = "test"))
  expect_error(proportion_S("test", also_single_tested = TRUE))

  # check too low amount of isolates
  expect_identical(
    suppressWarnings(proportion_R(example_isolates$AMX, minimum = nrow(example_isolates) + 1)),
    NA_real_
  )
  expect_identical(
    suppressWarnings(proportion_I(example_isolates$AMX, minimum = nrow(example_isolates) + 1)),
    NA_real_
  )
  expect_identical(
    suppressWarnings(proportion_S(example_isolates$AMX, minimum = nrow(example_isolates) + 1)),
    NA_real_
  )

  # warning for speed loss
  # expect_warning(proportion_R(as.character(example_isolates$GEN)))
  # expect_warning(proportion_I(as.character(example_isolates$GEN)))
  # expect_warning(proportion_S(example_isolates$AMC, as.character(example_isolates$GEN)))
  expect_error(proportion_df(c("A", "B", "C")))
  expect_error(proportion_df(example_isolates[, "date", drop = TRUE]))
})

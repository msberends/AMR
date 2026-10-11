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


test_that("test-sir_combine.R", {
  skip_on_cran()

  # truth table, written by hand: best result of two non-missing values
  vals <- c("S", "SDD", "I", "R", "WT", "NWT", "NS")
  best <- matrix(
    c(
      "S", "S", "S", "S", "S", "S", "S",
      "S", "SDD", "SDD", "SDD", "WT", "SDD", "SDD",
      "S", "SDD", "I", "I", "WT", "I", "I",
      "S", "SDD", "I", "R", "WT", "R", "R",
      "S", "WT", "WT", "WT", "WT", "WT", "WT",
      "S", "SDD", "I", "R", "WT", "NWT", "NWT",
      "S", "SDD", "I", "R", "WT", "NWT", "NS"
    ),
    nrow = 7, byrow = TRUE, dimnames = list(vals, vals)
  )
  expect_identical(best, t(best))
  # result of a non-missing value combined with a missing value (NA or NI), for only_all_tested = FALSE
  with_missing_eucast <- c(S = "S", SDD = "SDD", I = "I", R = NA, WT = "WT", NWT = NA, NS = NA)
  with_missing_clsi <- c(S = "S", SDD = "SDD", I = NA, R = NA, WT = "WT", NWT = NA, NS = NA)

  for (missing_val in c(NA, "NI")) {
    a <- c(rep(vals, each = length(vals)), vals, rep(missing_val, length(vals)), missing_val, "NI")
    b <- c(rep(vals, times = length(vals)), rep(missing_val, length(vals)), vals, NA, missing_val)
    expected_complete <- best[cbind(rep(vals, each = length(vals)), rep(vals, times = length(vals)))]
    expected_eucast <- c(expected_complete, with_missing_eucast[vals], with_missing_eucast[vals], NA, NA)
    expected_clsi <- c(expected_complete, with_missing_clsi[vals], with_missing_clsi[vals], NA, NA)
    expected_all_tested <- c(expected_complete, rep(NA, 2 * length(vals) + 2))

    a <- as.sir(a)
    b <- as.sir(b)
    expect_identical(sir_combine(a, b, guideline = "EUCAST"), as.sir(unname(expected_eucast)))
    expect_identical(sir_combine(a, b, guideline = "CLSI"), as.sir(unname(expected_clsi)))
    expect_identical(sir_combine(a, b, only_all_tested = TRUE, guideline = "EUCAST"), as.sir(unname(expected_all_tested)))
    expect_identical(sir_combine(a, b, only_all_tested = TRUE, guideline = "CLSI"), as.sir(unname(expected_all_tested)))
    # order of the antimicrobials does not matter
    expect_identical(sir_combine(a, b), sir_combine(b, a))
    expect_identical(sir_combine(a, b, guideline = "CLSI"), sir_combine(b, a, guideline = "CLSI"))
  }

  # three or more antimicrobials
  x1 <- as.sir(c("R", "R", "R", NA, "I", "R", "NI"))
  x2 <- as.sir(c("R", "I", "R", "R", "R", NA, "R"))
  x3 <- as.sir(c("R", "R", "S", "R", NA, "SDD", "R"))
  expect_identical(sir_combine(x1, x2, x3, guideline = "EUCAST"), as.sir(c("R", "I", "S", NA, "I", "SDD", NA)))
  expect_identical(sir_combine(x1, x2, x3, guideline = "CLSI"), as.sir(c("R", "I", "S", NA, NA, "SDD", NA)))
  expect_identical(sir_combine(x1, x2, x3, only_all_tested = TRUE), as.sir(c("R", "I", "S", NA, NA, NA, NA)))
  expect_identical(sir_combine(x1, x2, x3), sir_combine(x3, x1, x2))
  expect_identical(sir_combine(x1, x2, x3), sir_combine(sir_combine(x1, x2, only_all_tested = TRUE), x3))

  # input types
  expect_inherits(sir_combine(x1, x2), "sir")
  expect_identical(sir_combine(as.character(x1), as.character(x2)), sir_combine(x1, x2))
  expect_identical(sir_combine(data.frame(x1, x2, x3)), sir_combine(x1, x2, x3))
  expect_identical(sir_combine(data.frame(x1, x2), x3), sir_combine(x1, x2, x3))
  expect_identical(sir_combine(as.sir(character(0)), as.sir(character(0))), as.sir(character(0)))
  expect_identical(sir_combine(example_isolates$TZP, example_isolates$GEN), sir_combine(example_isolates[, c("TZP", "GEN")]))
  expect_identical(length(sir_combine(example_isolates$TZP, example_isolates$GEN)), nrow(example_isolates))
  expect_error(sir_combine(x1))
  expect_error(sir_combine(data.frame(x1)))
  expect_error(sir_combine())
  expect_error(sir_combine(x1, x2[1:3]))
  expect_error(sir_combine(x1, x2, only_all_tested = "yes"))
  expect_error(sir_combine(x1, x2, guideline = "something else"))

  # count and proportion functions use the same rules for combinations
  set.seed(321)
  n <- 500
  random_sir_values <- function() as.sir(sample(c("S", "SDD", "I", "R", "WT", "NWT", "NS", "NI", NA), n, replace = TRUE, prob = c(0.3, 0.05, 0.1, 0.25, 0.05, 0.05, 0.05, 0.05, 0.1)))
  r1 <- random_sir_values()
  r2 <- random_sir_values()
  r3 <- random_sir_values()
  for (oat in c(FALSE, TRUE)) {
    for (gl in c("EUCAST", "CLSI")) {
      combined <- sir_combine(r1, r2, r3, only_all_tested = oat, guideline = gl)
      expect_identical(count_all(r1, r2, r3, only_all_tested = oat), sum(!is.na(sir_combine(r1, r2, r3, only_all_tested = oat, guideline = "EUCAST"))))
      expect_identical(count_susceptible(r1, r2, r3, only_all_tested = oat, guideline = gl), count_susceptible(combined, guideline = gl))
      expect_identical(count_resistant(r1, r2, r3, only_all_tested = oat, guideline = gl), count_resistant(combined, guideline = gl))
      expect_equal(susceptibility(r1, r2, r3, only_all_tested = oat, guideline = gl, minimum = 0), susceptibility(combined, guideline = gl, minimum = 0))
      expect_equal(resistance(r1, r2, r3, only_all_tested = oat, guideline = gl, minimum = 0), resistance(combined, guideline = gl, minimum = 0))
      expect_equal(susceptibility(r1, r2, r3, only_all_tested = oat, guideline = gl, minimum = 0) + resistance(r1, r2, r3, only_all_tested = oat, guideline = gl, minimum = 0), 1)
      expect_identical(count_susceptible(combined, guideline = gl) + count_resistant(combined, guideline = gl), sum(!is.na(combined)))
    }
    combined <- sir_combine(r1, r2, only_all_tested = oat, guideline = "EUCAST")
    expect_identical(count_SI(r1, r2, only_all_tested = oat) + count_R(r1, r2, only_all_tested = oat), count_all(r1, r2, only_all_tested = oat))
    expect_identical(count_S(r1, r2, only_all_tested = oat), sum(combined %in% "S"))
    expect_identical(count_I(r1, r2, only_all_tested = oat), sum(combined %in% c("I", "SDD")))
    expect_identical(count_all(r1, r2, only_all_tested = oat), sum(!is.na(combined)))
    expect_identical(count_R(r1, r2, only_all_tested = oat), sum(combined %in% c("R", "NWT", "NS")))
    expect_equal(proportion_SI(r1, r2, only_all_tested = oat, minimum = 0) + proportion_R(r1, r2, only_all_tested = oat, minimum = 0), 1)
  }

  # antibiogram() gives the same result as the frozen implementation from before sir_combine() existed
  old_combination <- function(df, S_values, only_all_tested) {
    df[] <- lapply(df, function(col) {
      col <- as.sir(col)
      col[col == "NI"] <- NA_sir_
      col
    })
    x_transposed <- as.list(as.data.frame(t(df), stringsAsFactors = FALSE))
    if (isTRUE(only_all_tested)) {
      as.sir(vapply(FUN.VALUE = character(1), x_transposed, function(x) ifelse(anyNA(x), NA_character_, ifelse(any(x %in% S_values), "S", "R")), USE.NAMES = FALSE))
    } else {
      as.sir(vapply(FUN.VALUE = character(1), x_transposed, function(x) ifelse(any(x %in% S_values, na.rm = TRUE), "S", ifelse(anyNA(x), NA_character_, "R")), USE.NAMES = FALSE))
    }
  }
  df_random <- data.frame(mo = as.mo(sample(c("E. coli", "K. pneumoniae", "P. aeruginosa"), n, replace = TRUE)), TZP = r1, GEN = r2, TOB = r3, stringsAsFactors = FALSE)
  for (combine_SI in c(TRUE, FALSE)) {
    S_values <- c("S", "WT", if (combine_SI) c("SDD", "I"))
    for (oat in c(FALSE, TRUE)) {
      df_old <- df_random
      df_old$OLD <- old_combination(df_random[, c("TZP", "GEN", "TOB")], S_values, oat)
      ab_new <- suppressWarnings(suppressMessages(antibiogram(df_random, antimicrobials = "TZP+GEN+TOB", combine_SI = combine_SI, only_all_tested = oat, minimum = 0)))
      ab_old <- suppressWarnings(suppressMessages(antibiogram(df_old, antimicrobials = "OLD", combine_SI = combine_SI, only_all_tested = oat, minimum = 0)))
      long_new <- attributes(ab_new)$long_numeric
      long_old <- attributes(ab_old)$long_numeric
      expect_true(nrow(long_new) > 0)
      expect_identical(long_new[, colnames(long_new) != "ab"], long_old[, colnames(long_old) != "ab"])
      # wisca() uses the same combination step: the input of its Bayesian model must be identical
      wisca_new <- suppressWarnings(suppressMessages(wisca(df_random, antimicrobials = "TZP+GEN+TOB", combine_SI = combine_SI, only_all_tested = oat, simulations = 10, info = FALSE)))
      wisca_old <- suppressWarnings(suppressMessages(wisca(df_old, antimicrobials = "OLD", combine_SI = combine_SI, only_all_tested = oat, simulations = 10, info = FALSE)))
      params_new <- as.data.frame(attributes(wisca_new)$wisca_parameters)
      params_old <- as.data.frame(attributes(wisca_old)$wisca_parameters)
      expect_true(sum(params_new$n_tested) > 0)
      expect_identical(params_new[, colnames(params_new) != "ab"], params_old[, colnames(params_old) != "ab"])
    }
  }

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    out <- example_isolates %>% mutate(TZP_GEN = sir_combine(TZP, GEN))
    expect_identical(out$TZP_GEN, sir_combine(example_isolates$TZP, example_isolates$GEN))
    out <- example_isolates %>% mutate(any_aminoglycoside = sir_combine(GEN, TOB, AMK, KAN))
    expect_identical(out$any_aminoglycoside, sir_combine(example_isolates[, c("GEN", "TOB", "AMK", "KAN")]))
    expect_identical(
      example_isolates %>% group_by(ward) %>% summarise(s = susceptibility(TZP, GEN, minimum = 0)) %>% pull(s),
      example_isolates %>% mutate(TZP_GEN = sir_combine(TZP, GEN)) %>% group_by(ward) %>% summarise(s = susceptibility(TZP_GEN, minimum = 0)) %>% pull(s)
    )
  }
})

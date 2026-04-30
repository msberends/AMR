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

test_that("test-antibiogram.R", {
  skip_on_cran()

  # Traditional antibiogram ----------------------------------------------


  ab0 <- antibiogram(example_isolates)

  ab1 <- antibiogram(example_isolates,
    antimicrobials = c(aminoglycosides(), carbapenems())
  )

  ab2 <- antibiogram(example_isolates,
    antimicrobials = aminoglycosides(),
    ab_transform = "atc",
    mo_transform = "gramstain"
  )

  ab3 <- antibiogram(example_isolates,
    antimicrobials = carbapenems(),
    ab_transform = "ab",
    mo_transform = "name",
    formatting_type = 1
  )

  expect_inherits(ab1, "antibiogram")
  expect_inherits(ab2, "antibiogram")
  expect_inherits(ab3, "antibiogram")
  expect_equal(colnames(ab1), c("Pathogen", "Amikacin", "Gentamicin", "Imipenem", "Kanamycin", "Meropenem", "Tobramycin"))
  expect_equal(colnames(ab2), c("Pathogen", "J01GB01", "J01GB03", "J01GB04", "J01GB06"))
  expect_equal(colnames(ab3), c("Pathogen", "IPM", "MEM"))
  expect_equal(ab3$MEM, c(52, NA, 100, 100, NA))

  # Combined antibiogram -------------------------------------------------

  # combined antibiogram yield higher empiric coverage
  ab4 <- antibiogram(example_isolates,
    antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
    mo_transform = "gramstain"
  )

  ab5 <- antibiogram(example_isolates,
    antimicrobials = c("TZP", "TZP+TOB"),
    mo_transform = "gramstain",
    ab_transform = "name",
    sep = " & ",
    add_total_n = FALSE
  )

  expect_inherits(ab4, "antibiogram")
  expect_inherits(ab5, "antibiogram")
  expect_equal(colnames(ab4), c("Pathogen", "Piperacillin/tazobactam", "Piperacillin/tazobactam + Gentamicin", "Piperacillin/tazobactam + Tobramycin"))
  expect_equal(colnames(ab5), c("Pathogen", "Piperacillin/tazobactam", "Piperacillin/tazobactam & Tobramycin"))

  # Syndromic antibiogram ------------------------------------------------

  # the data set could contain a filter for e.g. respiratory specimens
  ab6 <- antibiogram(example_isolates,
    antimicrobials = c(aminoglycosides(), carbapenems()),
    syndromic_group = "ward",
    ab_transform = NULL
  )

  # with a custom language, though this will be determined automatically
  # (i.e., this table will be in Dutch on Dutch systems)
  ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
  ab7 <- antibiogram(ex1,
    antimicrobials = aminoglycosides(),
    ab_transform = "name",
    syndromic_group = ifelse(ex1$ward == "ICU",
      "IC", "Geen IC"
    ),
    language = "nl"
  )

  expect_inherits(ab6, "antibiogram")
  expect_inherits(ab7, "antibiogram")
  expect_equal(colnames(ab6), c("Syndromic Group", "Pathogen", "AMK", "GEN", "IPM", "KAN", "MEM", "TOB"))
  expect_equal(colnames(ab7), c("Syndroomgroep", "Pathogeen", "Amikacine", "Gentamicine", "Tobramycine"))

  # Weighted-incidence syndromic combination antibiogram (WISCA) ---------

  # the data set could contain a filter for e.g. respiratory specimens
  ab8 <- suppressWarnings(antibiogram(example_isolates,
    antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
    wisca = TRUE
  ))

  expect_inherits(ab8, "antibiogram")
  expect_inherits(retrieve_wisca_parameters(ab8), "data.frame")
  expect_inherits(attributes(ab8)$long_numeric, "data.frame")
  expect_equal(colnames(ab8), c("Piperacillin/tazobactam", "Piperacillin/tazobactam + Gentamicin", "Piperacillin/tazobactam + Tobramycin"))

  # grouped tibbles
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_warning(
      ab9 <- example_isolates %>%
        group_by(ward, gender) %>%
        wisca(antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
    )
    expect_equal(colnames(ab9), c("ward", "gender", "Piperacillin/tazobactam", "Piperacillin/tazobactam + Gentamicin", "Piperacillin/tazobactam + Tobramycin"))
  }

  # Parallel computing ----------------------------------------------------
  # Tests must pass even when only 1 core is available; parallel = TRUE then
  # silently falls back to sequential, but results must still be correct.

  if (AMR:::pkg_is_available("future.apply")) {
    set.seed(42)

    # sequential reference for WISCA
    wisca_seq <- suppressWarnings(suppressMessages(
      wisca(example_isolates, antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"), simulations = 100, info = FALSE)
    ))

    future::plan(future::multicore)

    # 1. parallel = TRUE produces the same antibiogram structure as sequential
    wisca_par <- suppressWarnings(suppressMessages(
      wisca(example_isolates, antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"), simulations = 100, parallel = TRUE, info = FALSE)
    ))
    expect_inherits(wisca_par, "antibiogram")
    expect_equal(colnames(wisca_par), colnames(wisca_seq))
    expect_true(isTRUE(attributes(wisca_par)$wisca))

    # 2. coverage values fall within [0, 100] (basic sanity)
    ln <- attributes(wisca_par)$long_numeric
    expect_true(all(ln$coverage >= 0 & ln$coverage <= 1, na.rm = TRUE))
    expect_true(all(ln$lower_ci <= ln$coverage, na.rm = TRUE))
    expect_true(all(ln$upper_ci >= ln$coverage, na.rm = TRUE))

    # 3. a second parallel run gives the same column names
    wisca_par2 <- suppressWarnings(suppressMessages(
      wisca(example_isolates, antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"), simulations = 100, parallel = TRUE, info = FALSE)
    ))
    expect_equal(colnames(wisca_par), colnames(wisca_par2))

    # 4. parallel with workers = 1 gives same structure as sequential
    future::plan(future::multicore, workers = 1)
    wisca_par1 <- suppressWarnings(suppressMessages(
      wisca(example_isolates, antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"), simulations = 100, parallel = TRUE, info = FALSE)
    ))
    expect_equal(colnames(wisca_seq), colnames(wisca_par1))

    # 5. grouped antibiogram in parallel yields identical structure to sequential
    if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
      future::plan(future::sequential)
      ab_grp_seq <- suppressWarnings(suppressMessages(
        example_isolates %>%
          group_by(ward) %>%
          wisca(antimicrobials = c("TZP", "TZP+TOB"), simulations = 50, info = FALSE)
      ))
      future::plan(future::multicore)
      ab_grp_par <- suppressWarnings(suppressMessages(
        example_isolates %>%
          group_by(ward) %>%
          wisca(antimicrobials = c("TZP", "TZP+TOB"), simulations = 50, parallel = TRUE, info = FALSE)
      ))
      expect_equal(colnames(ab_grp_seq), colnames(ab_grp_par))
      expect_equal(nrow(ab_grp_seq), nrow(ab_grp_par))
    }

    # 6. parallel = TRUE without a plan raises an informative error
    future::plan(future::sequential)
    expect_error(
      suppressWarnings(wisca(example_isolates, antimicrobials = "TZP", parallel = TRUE, info = FALSE)),
      "non-sequential"
    )

    future::plan(future::sequential)
  }

  # Generate plots with ggplot2 or base R --------------------------------

  pdf(NULL) # prevent Rplots.pdf being created

  expect_silent(plot(ab1))
  expect_silent(plot(ab2))
  expect_silent(plot(ab3))
  expect_silent(plot(ab4))
  expect_silent(plot(ab5))
  expect_silent(plot(ab6))
  expect_silent(plot(ab7))
  expect_silent(plot(ab8))
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_silent(plot(ab9))
  }

  if (AMR:::pkg_is_available("ggplot2")) {
    expect_inherits(ggplot2::autoplot(ab1), "gg")
    expect_inherits(ggplot2::autoplot(ab2), "gg")
    expect_inherits(ggplot2::autoplot(ab3), "gg")
    expect_inherits(ggplot2::autoplot(ab4), "gg")
    expect_inherits(ggplot2::autoplot(ab5), "gg")
    expect_inherits(ggplot2::autoplot(ab6), "gg")
    expect_inherits(ggplot2::autoplot(ab7), "gg")
    expect_inherits(ggplot2::autoplot(ab8), "gg")
    if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
      expect_inherits(ggplot2::autoplot(ab9), "gg")
    }
  }
})

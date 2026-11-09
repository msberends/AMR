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

test_that("test-ggplot_sir.R", {
  skip_on_cran()

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE) &&
    AMR:::pkg_is_available("ggplot2", also_load = TRUE)) {
    pdf(NULL) # prevent Rplots.pdf being created

    # data should be equal
    expect_equal(
      (example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir())$data %>%
        summarise_all(resistance) %>%
        as.double(),
      example_isolates %>%
        select(AMC, CIP) %>%
        summarise_all(resistance) %>%
        as.double()
    )

    expect_inherits(
      example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "interpretation", facet = "antibiotic"),
      "gg"
    )
    expect_inherits(
      example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "antibiotic", facet = "interpretation"),
      "gg"
    )
    # the facet is passed on as value, so the plot can be built
    expect_inherits(
      ggplot_build(example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "interpretation", facet = "antibiotic")),
      "ggplot_built"
    )

    expect_equal(
      (example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "interpretation", facet = "antibiotic"))$data %>%
        summarise_all(resistance) %>%
        as.double(),
      example_isolates %>%
        select(AMC, CIP) %>%
        summarise_all(resistance) %>%
        as.double()
    )

    expect_equal(
      (example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "antibiotic", facet = "interpretation"))$data %>%
        summarise_all(resistance) %>%
        as.double(),
      example_isolates %>%
        select(AMC, CIP) %>%
        summarise_all(resistance) %>%
        as.double()
    )

    expect_equal(
      (example_isolates %>%
        select(AMC, CIP) %>%
        ggplot_sir(x = "antibiotic", facet = "interpretation"))$data %>%
        summarise_all(count_resistant) %>%
        as.double(),
      example_isolates %>%
        select(AMC, CIP) %>%
        summarise_all(count_resistant) %>%
        as.double()
    )

    # the 'I' category is merged with R for CLSI
    expect_identical(
      levels(ggplot_sir(example_isolates[, c("AMC", "CIP")], guideline = "CLSI")$layers[[1]]$data(example_isolates[, c("AMC", "CIP")])$interpretation),
      c("S", "IR")
    )
    # S, I, and R are kept separate with guideline = NULL
    expect_identical(
      levels(ggplot_sir(example_isolates[, c("AMC", "CIP")], guideline = NULL)$layers[[1]]$data(example_isolates[, c("AMC", "CIP")])$interpretation),
      c("S", "I", "R")
    )
    # the deprecated combine_SI must be followed literally if set, regardless of the guideline
    expect_identical(
      levels(suppressWarnings(ggplot_sir(example_isolates[, c("AMC", "CIP")], combine_SI = TRUE, guideline = "CLSI"))$layers[[1]]$data(example_isolates[, c("AMC", "CIP")])$interpretation),
      c("SI", "R")
    )

    # support for scale_type ab and mo
    expect_inherits(
      (data.frame(
        mo = as.mo(c("e. coli", "s aureus")),
        n = c(40, 100)
      ) %>%
        ggplot(aes(x = mo, y = n)) +
        geom_col())$data,
      "data.frame"
    )
    expect_inherits(
      (data.frame(
        ab = as.ab(c("amx", "amc")),
        n = c(40, 100)
      ) %>%
        ggplot(aes(x = ab, y = n)) +
        geom_col())$data,
      "data.frame"
    )

    expect_inherits(
      (data.frame(
        ab = as.ab(c("amx", "amc")),
        n = c(40, 100)
      ) %>%
        ggplot(aes(x = ab, y = n)) +
        geom_col())$data,
      "data.frame"
    )

    # support for manual colours
    expect_inherits(
      suppressWarnings((ggplot(data.frame(
        x = c("Value1", "Value2", "Value3"),
        y = c(1, 2, 3),
        z = c("Value4", "Value5", "Value6")
      )) +
        geom_col(aes(x = x, y = y, fill = z)) +
        scale_sir_colours(aesthetics = "fill", Value4 = "S", Value5 = "I", Value6 = "R"))$data),
      "data.frame"
    )

    # datalabels are placed on their own bar, also when dodged or stacked
    for (pos in c("fill", "stack", "dodge")) {
      built <- ggplot_build(ggplot_sir(example_isolates[, c("AMX", "CIP")], position = pos, guideline = NULL))
      bars <- built$data[[1]]
      labels <- built$data[[2]]
      bars <- bars[order(bars$x, bars$group), , drop = FALSE]
      labels <- labels[order(labels$x, labels$group), , drop = FALSE]
      expect_equal(labels$x, bars$x)
      if (pos == "dodge") {
        expect_equal(labels$y, bars$ymax)
      } else {
        expect_equal(labels$y, (bars$ymin + bars$ymax) / 2)
      }
    }

    # without the interpretation in the plot, only the susceptible proportion is shown, following the guideline
    grouped <- example_isolates[, c("gender", "AMX")]
    grouped <- structure(grouped, class = c("grouped_df", "tbl_df", "tbl", "data.frame"), groups = data.frame(gender = c("F", "M"), .rows = I(list(which(grouped$gender == "F"), which(grouped$gender == "M")))))
    for (gl in list("EUCAST", "CLSI", NULL)) {
      expected <- sir_df(grouped, guideline = gl)
      expected <- expected[which(as.character(expected$interpretation) == levels(expected$interpretation)[1]), , drop = FALSE]
      p <- ggplot_sir(grouped, x = "antibiotic", fill = "gender", position = "dodge", guideline = gl)
      built <- ggplot_build(p)
      expect_equal(sort(built$data[[1]]$y), sort(expected$value))
      expect_equal(nrow(built$data[[2]]), nrow(expected))
      # data labels show the number of tested isolates
      expected_n <- tapply(sir_df(grouped, guideline = gl)$isolates, sir_df(grouped, guideline = gl)$gender, sum)
      expect_identical(sort(built$data[[2]]$label), sort(paste0("n=", expected_n)))
      expect_identical(p$labels$y, paste("Proportion", levels(expected$interpretation)[1]))
    }
    # all interpretations when they are shown in the facets
    built <- ggplot_build(ggplot_sir(grouped, x = "antibiotic", fill = "gender", facet = "interpretation", position = "dodge", guideline = "EUCAST"))
    expect_equal(nrow(built$data[[1]]), nrow(sir_df(grouped, guideline = "EUCAST")))

    # split by antibiotic if it is not shown in the plot and there are several
    grouped2 <- example_isolates[, c("gender", "AMX", "CIP")]
    grouped2 <- structure(grouped2, class = c("grouped_df", "tbl_df", "tbl", "data.frame"), groups = data.frame(gender = c("F", "M"), .rows = I(list(which(grouped2$gender == "F"), which(grouped2$gender == "M")))))
    p <- ggplot_sir(grouped2, x = "gender", guideline = "EUCAST")
    expect_true(inherits(p$facet, "FacetWrap"))
    expect_identical(names(p$facet$params$facets), "antibiotic")
    p <- ggplot_sir(grouped2, x = "interpretation", fill = "interpretation", facet = "gender", guideline = "EUCAST")
    expect_true(inherits(p$facet, "FacetGrid"))
    expect_identical(names(p$facet$params$cols), "antibiotic")
    expect_identical(names(p$facet$params$rows), "gender")
    # not if the antibiotic is shown, or if there is only one
    expect_false(inherits(ggplot_sir(grouped2, x = "antibiotic", guideline = "EUCAST")$facet, "FacetWrap"))
    expect_false(inherits(ggplot_sir(grouped, x = "gender", guideline = "EUCAST")$facet, "FacetWrap"))
  }
})

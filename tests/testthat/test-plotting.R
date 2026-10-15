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

test_that("test-plotting.R", {
  skip_on_cran()

  if (AMR:::pkg_is_available("ggplot2", also_load = TRUE)) {
    pdf(NULL) # prevent Rplots.pdf being created

    # scale_*_mic
    aesthetics <- c("x", "y", "colour", "fill")
    expected_methods <- c("transform", "transform_df", "breaks", "labels", "limits")
    for (aest in aesthetics) {
      scale_fn_name <- paste0("scale_", aest, "_continuous")
      scale_obj <- getExportedValue("ggplot2", scale_fn_name)()
      for (method in expected_methods) {
        expect_true(is.function(scale_obj[[method]]) || method %in% names(scale_obj),
          info = paste0("Method '", method, "' is missing in ggplot2::", scale_fn_name)
        )
      }
    }

    # scale_*_sir
    aesthetics <- c("colour", "fill")
    expected_methods <- c("transform", "transform_df", "labels", "limits")
    for (aest in aesthetics) {
      scale_fn_name <- "scale_discrete_manual"
      scale_obj <- getExportedValue("ggplot2", scale_fn_name)(aesthetics = aest)
      for (method in expected_methods) {
        expect_true(is.function(scale_obj[[method]]) || method %in% names(scale_obj),
          info = paste0("Method '", method, "' is missing in ggplot2::", scale_fn_name)
        )
      }
    }
    for (method in expected_methods) {
      expect_true(is.function(ggplot2::scale_x_discrete()[[method]]) || method %in% names(ggplot2::scale_x_discrete()),
        info = paste0("Method '", method, "' is missing in ggplot2::", "scale_x_discrete")
      )
    }

    expect_inherits(
      ggplot(
        data.frame(
          count = c(1, 2, 3, 4),
          sir = c("S", "I", "R", "SDD")
        ),
        aes(x = sir, y = count, fill = sir)
      ) +
        geom_col() +
        scale_x_sir(eucast_I = F, language = "el") +
        scale_fill_sir(eucast_I = T, language = "nl"),
      "gg"
    )
    expect_inherits(
      ggplot(
        data.frame(
          mic = as.mic(c(2, 4, 8, 16)),
          sir = as.sir(c("S", "I", "R", "SDD"))
        ),
        aes(x = sir, y = mic)
      ) +
        geom_point() +
        scale_y_mic(),
      "gg"
    )
    expect_inherits(
      ggplot(
        data.frame(
          mic = as.mic(c(2, 4, 8, 16)),
          sir = as.sir(c("S", "I", "R", "SDD"))
        ),
        aes(x = sir, y = mic)
      ) +
        geom_col() +
        scale_y_mic(),
      "gg"
    )
    expect_inherits(
      ggplot(
        data.frame(
          mic = as.mic(c(2, 4, 8, 16)),
          sir = as.sir(c("S", "I", "R", "SDD"))
        ),
        aes(x = sir, y = mic)
      ) +
        geom_col() +
        scale_y_mic(mic_range = c(4, 16)) +
        scale_x_sir(),
      "gg"
    )

    # named colours are matched by name, and SDD, S/I and I/R follow S and R
    expect_identical(AMR:::expand_SIR_colours(c(R = "r", I = "i", SDD = "d", S = "s")), c("s", "d", "i", "r", "s", "r"))
    expect_identical(AMR:::expand_SIR_colours(c(S = "s", I = "i", R = "r")), c("s", "s", "i", "r", "s", "r"))
    expect_identical(AMR:::expand_SIR_colours(c("s", "i", "r")), c("s", "s", "i", "r", "s", "r"))
    expect_error(AMR:::expand_SIR_colours(c(S = "s", I = "i", X = "r")))
    # missing names S, I and R get the default colours, as in the example of ggplot_sir()
    expect_identical(AMR:::expand_SIR_colours(c(SI = "y")), c("#3CAEA3", "#3CAEA3", "#F6D55C", "#ED553B", "y", "#ED553B"))

    # colours of SIR bar plots follow the SIR values
    sir_values <- as.sir(c("S", "S", "R", "NI", "NS"))
    keep <- AMR:::plot_sir_levels_to_keep(sir_values, table(sir_values))
    expect_identical(names(table(sir_values))[keep], c("S", "R", "NI", "NS"))
    p <- autoplot(sir_values)
    expect_identical(as.character(p$data$x), c("S", "R", "NI", "NS"))
    expect_identical(unique(ggplot_build(p)$data[[1]]$fill), c("#3CAEA3", "#ED553B", "grey30"))

    # SDD in disk diffusion plots
    p <- autoplot(as.disk(15:30), mo = "Escherichia coli", ab = "cefepime", guideline = "CLSI 2025", language = "en")
    expect_true("(SDD) Susceptible dose-dependent" %in% as.character(p$data$cols))
    expect_false(anyNA(p$data$cols))

    # colours and labels are named by, and cover, all SIR levels
    expect_identical(names(AMR:::plot_sir_colour_values(AMR:::expand_SIR_colours("red"))), AMR:::VALID_SIR_LEVELS)
    expect_identical(names(AMR:::plot_sir_label_values("EUCAST")), AMR:::VALID_SIR_LEVELS)

    # arguments are passed on
    expect_identical(scale_x_sir(name = "Interpretation")$name, "Interpretation")
    expect_identical(scale_fill_sir(values = c(S = "black"))$palette(1)[["S"]], "black")
    expect_identical(scale_fill_sir(values = c(S = "black"))$palette(1)[["R"]], "#ED553B")
    expect_error(scale_fill_sir(values = c("black", "white", "grey")))
    expect_warning(scale_x_sir(values = c(S = "black")))
    expect_identical(AMR:::plot_count_ticks(3), c(0, 1, 2, 3))
    expect_identical(AMR:::plot_name_of_I("EUCAST"), "Susceptible, incr. exp.")
    expect_identical(AMR:::plot_name_of_I("EUCAST 2018"), "Intermediate")

    # unquoted column names in facet_sir() are not evaluated, also not if an object with that name exists
    antibiotic <- "some other column"
    expect_identical(facet_sir(antibiotic)$params$facets, facet_sir("antibiotic")$params$facets)

    # the default of eucast_I follows the guideline option
    old_option <- options(AMR_guideline = "EUCAST")
    expect_identical(scale_fill_sir(language = "en")$labels("I"), "(I) Susceptible, incr. exp.")
    options(AMR_guideline = "EUCAST 2018")
    expect_identical(scale_fill_sir(language = "en")$labels("I"), "(I) Intermediate")
    options(AMR_guideline = "CLSI")
    expect_identical(scale_fill_sir(language = "en")$labels("I"), "(I) Intermediate")
    options(old_option)

    # keep_operators = TRUE/FALSE equals "all"/"none"
    mic_labels <- function(keep_operators) {
      p <- ggplot(data.frame(m = as.mic(c(">=8", "2", "<=0.5"))), aes(m, 1)) +
        geom_point() +
        scale_x_mic(keep_operators = keep_operators)
      ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
    }
    expect_identical(mic_labels(TRUE), mic_labels("all"))
    expect_identical(mic_labels(FALSE), mic_labels("none"))
    expect_false(identical(mic_labels(TRUE), mic_labels(FALSE)))
  }

  # base R plots: texts in the legend (drawn with text()) and bar names (drawn with axis())
  recorded_strings <- function(expr, fn) {
    pdf(NULL)
    dev.control("enable")
    on.exit(dev.off())
    force(expr)
    unlist(lapply(recordPlot()[[1]], function(call) {
      if (identical(call[[2]][[1]]$name, fn)) {
        as.character(call[[2]][[if (fn == "C_text") 3 else 4]])
      }
    }))
  }
  mic_values <- as.mic(c(0.25, 0.25, 4, 32))

  # the legend is based on the interpretation, so SDD is not shown when it has the colour of S
  expect_identical(
    recorded_strings(plot(mic_values, mo = "Escherichia coli", ab = "AMC", language = "en", colours_SIR = c("green", "yellow", "red")), "C_text"),
    c("(S) Susceptible", "(R) Resistant")
  )
  # zone diameters run from resistant to susceptible, and so does the legend
  expect_identical(
    recorded_strings(plot(as.disk(15:30), mo = "Escherichia coli", ab = "cefepime", guideline = "CLSI 2025", language = "en"), "C_text"),
    c("(R) Resistant", "(SDD) Susceptible dose-dependent", "(S) Susceptible")
  )
  # barplot() for MIC values passes on language and expand
  expect_identical(
    recorded_strings(barplot(mic_values, mo = "Escherichia coli", ab = "AMC", language = "nl"), "C_text"),
    c("(S) Gevoelig", "(R) Resistent")
  )
  expect_true(length(recorded_strings(barplot(mic_values, expand = TRUE), "C_axis")) > length(recorded_strings(barplot(mic_values, expand = FALSE), "C_axis")))
})

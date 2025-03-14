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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

test_that("test-plotting.R", {
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
  }
})

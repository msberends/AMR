# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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

sir <- random_sir(100)
rsi <- sir
class(rsi) <- gsub("sir", "rsi", class(rsi))
mic <- random_mic(100)
disk <- random_disk(100)

expect_identical(summary(sir), summary(rsi))
expect_identical(c(sir), c(rsi))

expect_identical(suppressWarnings(suppressMessages(as.rsi(as.character(rsi)))),
                 suppressWarnings(suppressMessages(as.sir(as.character(sir)))))
expect_identical(suppressWarnings(suppressMessages(as.rsi(mic, mo = "Escherichia coli", ab = "CIP"))),
                 suppressWarnings(suppressMessages(as.sir(mic, mo = "Escherichia coli", ab = "CIP"))))
expect_identical(suppressWarnings(suppressMessages(as.rsi(disk, mo = "Escherichia coli", ab = "CIP"))),
                 suppressWarnings(suppressMessages(as.sir(disk, mo = "Escherichia coli", ab = "CIP"))))
expect_identical(suppressWarnings(suppressMessages(as.rsi(data.frame(CIP = mic, mo = "Escherichia coli")))),
                 suppressWarnings(suppressMessages(as.sir(data.frame(CIP = mic, mo = "Escherichia coli")))))

expect_identical(suppressWarnings(n_rsi(example_isolates$CIP)),
                 suppressWarnings(n_sir(example_isolates$CIP)))

expect_identical(suppressWarnings(rsi_df(example_isolates)),
                 suppressWarnings(sir_df(example_isolates)))

expect_identical(suppressWarnings(is.rsi.eligible(example_isolates)),
                 suppressWarnings(is_sir_eligible(example_isolates)))

if (AMR:::pkg_is_available("ggplot2")) {
  expect_equal(suppressWarnings(ggplot_rsi(example_isolates[, c("CIP", "GEN", "TOB")])),
               suppressWarnings(ggplot_sir(example_isolates[, c("CIP", "GEN", "TOB")])))
  
  p <- ggplot2::ggplot(example_isolates[, c("CIP", "GEN", "TOB")])
  expect_equal(suppressWarnings(p + geom_rsi() + scale_rsi_colours() + labels_rsi_count() + facet_rsi() + theme_rsi()),
               suppressWarnings(p + geom_sir() + scale_sir_colours() + labels_sir_count() + facet_sir() + theme_sir()))
}

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

  expect_stdout(print(example_isolates %>%
    select(AMC, CIP) %>%
    ggplot_sir(x = "interpretation", facet = "antibiotic")))
  expect_stdout(print(example_isolates %>%
    select(AMC, CIP) %>%
    ggplot_sir(x = "antibiotic", facet = "interpretation")))

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
    (ggplot(data.frame(
      x = c("Value1", "Value2", "Value3"),
      y = c(1, 2, 3),
      z = c("Value4", "Value5", "Value6")
    )) +
      geom_col(aes(x = x, y = y, fill = z)) +
      scale_sir_colours(Value4 = "S", Value5 = "I", Value6 = "R"))$data,
    "data.frame"
  )
}

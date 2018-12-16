# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

context("ggplot_rsi.R")

test_that("ggplot_rsi works", {

  skip_if_not("ggplot2" %in% rownames(installed.packages()))

  library(dplyr)
  library(ggplot2)

  # data should be equal
  expect_equal(
    (septic_patients %>% select(amcl, cipr) %>% ggplot_rsi())$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(amcl, cipr) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(amcl, cipr) %>% ggplot_rsi(x = "Interpretation", facet = "Antibiotic"))$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(amcl, cipr) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(amcl, cipr) %>% ggplot_rsi(x = "Antibiotic", facet = "Interpretation"))$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(amcl, cipr) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(amcl, cipr) %>% ggplot_rsi(x = "Antibiotic",
                                                           facet = "Interpretation",
                                                           fun = count_df))$data %>%
      summarise_all(count_IR) %>% as.double(),
    septic_patients %>% select(amcl, cipr) %>%
      summarise_all(count_IR) %>% as.double()
  )

  expect_equal(colnames(getlbls(septic_patients %>% select(amcl, cipr))),
               c("Interpretation", "Antibiotic", "Value", "lbl"))

  expect_error(ggplot_rsi(septic_patients, fun = "invalid"))
  expect_error(geom_rsi(septic_patients, fun = "invalid"))

})

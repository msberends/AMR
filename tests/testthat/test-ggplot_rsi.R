# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("ggplot_rsi.R")

test_that("ggplot_rsi works", {

  skip_if_not("ggplot2" %in% rownames(installed.packages()))

  library(dplyr)
  library(ggplot2)

  # data should be equal
  expect_equal(
    (septic_patients %>% select(AMC, CIP) %>% ggplot_rsi())$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(AMC, CIP) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(AMC, CIP) %>% ggplot_rsi(x = "Interpretation", facet = "Antibiotic"))$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(AMC, CIP) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(AMC, CIP) %>% ggplot_rsi(x = "Antibiotic", facet = "Interpretation"))$data %>%
      summarise_all(portion_IR) %>% as.double(),
    septic_patients %>% select(AMC, CIP) %>%
      summarise_all(portion_IR) %>% as.double()
  )

  expect_equal(
    (septic_patients %>% select(AMC, CIP) %>% ggplot_rsi(x = "Antibiotic",
                                                           facet = "Interpretation",
                                                           fun = count_df))$data %>%
      summarise_all(count_IR) %>% as.double(),
    septic_patients %>% select(AMC, CIP) %>%
      summarise_all(count_IR) %>% as.double()
  )

  expect_error(ggplot_rsi(septic_patients, fun = "invalid"))
  expect_error(geom_rsi(septic_patients, fun = "invalid"))

  # support for scale_type ab and mo
  expect_equal(class((data.frame(mo = as.mo(c("e. coli", "s aureus")),
                                 n = c(40, 100)) %>%
                        ggplot(aes(x = mo, y = n)) +
                        geom_col())$data),
               "data.frame")
  expect_equal(class((data.frame(ab = as.ab(c("amx", "amc")),
                                 n = c(40, 100)) %>%
                        ggplot(aes(x = ab, y = n)) +
                        geom_col())$data),
               "data.frame")

})

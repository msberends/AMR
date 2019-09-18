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

context("mdro.R")

test_that("mdro works", {
  library(dplyr)

  expect_error(suppressWarnings(mdro(example_isolates, country = "invalid", col_mo = "mo", info = TRUE)))
  expect_error(suppressWarnings(mdro(example_isolates, country = "fr", info = TRUE)))
  expect_error(mdro(example_isolates, guideline = c("BRMO", "MRGN"), info = TRUE))
  expect_error(mdro(example_isolates, col_mo = "invalid", info = TRUE))

  outcome <- mdro(example_isolates)
  outcome <- eucast_exceptional_phenotypes(example_isolates, info = TRUE)
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  outcome <- mdro(example_isolates, "nl", info = TRUE)
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  # example_isolates should have these finding using Dutch guidelines
  expect_equal(outcome %>% freq() %>% pull(count),
               c(1972, 22, 6)) # 1969 neg, 25 unconfirmed, 6 pos

  expect_equal(brmo(example_isolates, info = FALSE),
               mdro(example_isolates, guideline = "BRMO", info = FALSE))

  # still working on German guidelines
  expect_error(suppressWarnings(mrgn(example_isolates, info = TRUE)))

  # test Dutch P. aeruginosa MDRO
  expect_equal(
    as.character(mdro(data.frame(mo = as.mo("P. aeruginosa"),
                                 cfta = "S",
                                 cipr = "S",
                                 mero = "S",
                                 imip = "S",
                                 gent = "S",
                                 tobr = "S",
                                 pita = "S"),
                      guideline = "BRMO",
                      col_mo = "mo",
                      info = FALSE)),
    "Negative")
  expect_equal(
    as.character(mdro(data.frame(mo = as.mo("P. aeruginosa"),
                                 cefta = "R",
                                 cipr = "R",
                                 mero = "R",
                                 imip = "R",
                                 gent = "R",
                                 tobr = "R",
                                 pita = "R"),
                      guideline = "BRMO",
                      col_mo = "mo",
                      info = FALSE)),
    "Positive")

  # MDR TB
  expect_equal(
    # select only rifampicine, mo will be determined automatically (as M. tuberculosis),
    # number of mono-resistant strains should be equal to number of rifampicine-resistant strains
    example_isolates %>% select(RIF) %>% mdr_tb() %>% freq() %>% pull(count) %>% .[2],
    count_R(example_isolates$RIF))

  sample_rsi <- function() {
    sample(c("S", "I", "R"),
           size = 5000,
           prob = c(0.5, 0.1, 0.4),
           replace = TRUE)
  }
  expect_gt(
    #suppressWarnings(
      data.frame(rifampicin = sample_rsi(),
                 inh = sample_rsi(),
                 gatifloxacin = sample_rsi(),
                 eth = sample_rsi(),
                 pza = sample_rsi(),
                 MFX = sample_rsi(),
                 KAN = sample_rsi()) %>%
        mdr_tb() %>%
        n_distinct()
      #)
      ,
    2)

})

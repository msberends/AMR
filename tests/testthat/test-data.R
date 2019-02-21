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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

context("data.R")

test_that("data sets are valid", {
  # IDs should always be unique
  expect_identical(nrow(antibiotics), length(unique(antibiotics$atc)))
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))

  # there should be no diacritics (i.e. non ASCII) characters in the datasets
  library(dplyr)
  # check only character variables:
  test_microorganisms <- microorganisms %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_microorganisms.old <- microorganisms.old %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_antibiotics <- antibiotics %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_septic_patients <- septic_patients %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  # and compare them with their transformed version:
  expect_identical(test_microorganisms,
                   test_microorganisms %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_microorganisms.old,
                   test_microorganisms.old %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_antibiotics,
                   test_antibiotics %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_septic_patients,
                   test_septic_patients %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

})


test_that("creation of data sets is valid", {
  DT <- make_DT()
  expect_lt(nrow(DT[prevalence == 1]), nrow(DT[prevalence == 2]))
  expect_lt(nrow(DT[prevalence == 2]), nrow(DT[prevalence == 3]))
})

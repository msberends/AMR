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

context("data.R")

test_that("data sets are valid", {
  # IDs should always be unique
  expect_identical(nrow(antibiotics), length(unique(antibiotics$atc)))
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))

  # there should be no diacritics (i.e. non ASCII) characters in the datasets
  datasets <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item"]
  for (i in 1:length(datasets)) {
    dataset <- get(datasets[i], envir = asNamespace("AMR"))
    expect_identical(dataset_UTF8_to_ASCII(dataset), dataset)
  }
})

test_that("creation of data sets is valid", {
  df <- make()
  expect_lt(nrow(df[which(df$prevalence == 1), ]), nrow(df[which(df$prevalence == 2), ]))
  expect_lt(nrow(df[which(df$prevalence == 2), ]), nrow(df[which(df$prevalence == 3), ]))
  DT <- make_DT()
  expect_lt(nrow(DT[prevalence == 1]), nrow(DT[prevalence == 2]))
  expect_lt(nrow(DT[prevalence == 2]), nrow(DT[prevalence == 3]))
  old <- make_trans_tbl()
  expect_gt(length(old), 0)
})

test_that("CoL version info works", {
 expect_identical(class(catalogue_of_life_version()),
                  c("catalogue_of_life_version", "list"))
})

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("data.R")

test_that("data sets are valid", {
  
  expect_true(check_dataset_integrity()) # in misc.R
  
  # IDs should always be unique
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))
  expect_identical(class(microorganisms$mo), "mo")
  expect_identical(nrow(antibiotics), length(unique(antibiotics$ab)))
  expect_identical(class(antibiotics$ab), "ab")
  
  # check cross table reference
  expect_true(all(microorganisms.codes$mo %in% microorganisms$mo))
  expect_true(all(example_isolates$mo %in% microorganisms$mo))
  expect_true(all(microorganisms.translation$mo_new %in% microorganisms$mo))
  expect_true(all(rsi_translation$mo %in% microorganisms$mo))
  expect_false(any(is.na(microorganisms.codes$code)))
  expect_false(any(is.na(microorganisms.codes$mo)))

  # antibiotic names must always be coercible to their original AB code
  expect_identical(antibiotics$ab, as.ab(antibiotics$name))

  # there should be no diacritics (i.e. non ASCII) characters in the datasets (CRAN policy)
  datasets <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item"]
  for (i in seq_len(length(datasets))) {
    dataset <- get(datasets[i], envir = asNamespace("AMR"))
    expect_identical(dataset_UTF8_to_ASCII(dataset), dataset, label = datasets[i])
  }
})

test_that("creation of data sets is valid", {
  DT <- make_DT()
  expect_lt(nrow(DT[prevalence == 1]), nrow(DT[prevalence == 2]))
  expect_lt(nrow(DT[prevalence == 2]), nrow(DT[prevalence == 3]))
  expect_true(all(c("mo", "fullname",
                    "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies",
                    "rank", "col_id", "species_id", "source", "ref", "prevalence",
                    "kingdom_index", "fullname_lower", "g_species") %in% colnames(DT)))

  oldDT <- make_oldDT()
  expect_true(all(c("col_id", "col_id_new", "fullname", "ref", "prevalence",
                    "fullname_lower", "g_species") %in% colnames(oldDT)))
  
  old <- make_trans_tbl()
  expect_gt(length(old), 0)
  
})

test_that("CoL version info works", {
 expect_identical(class(catalogue_of_life_version()),
                  c("catalogue_of_life_version", "list"))

  expect_output(print(catalogue_of_life_version()))
})

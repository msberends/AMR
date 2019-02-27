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

context("eucast_rules.R")

test_that("EUCAST rules work", {

  expect_error(suppressWarnings(eucast_rules(septic_patients, col_mo = "Non-existing")))

  expect_identical(colnames(septic_patients),
                   colnames(suppressWarnings(eucast_rules(septic_patients))))

  a <- data.frame(mo = c("Klebsiella pneumoniae",
                         "Pseudomonas aeruginosa",
                         "Enterobacter aerogenes"),
                  amox = "-",           # Amoxicillin
                  stringsAsFactors = FALSE)
  b <- data.frame(mo = c("Klebsiella pneumoniae",
                         "Pseudomonas aeruginosa",
                         "Enterobacter aerogenes"),
                  amox = "R",       # Amoxicillin
                  stringsAsFactors = FALSE)
  expect_identical(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)
  expect_identical(suppressWarnings(eucast_rules(a, "mo", info = TRUE)), b)
  expect_identical(suppressWarnings(interpretive_reading(a, "mo", info = TRUE)), b)

  a <- data.frame(mo = c("Staphylococcus aureus",
                         "Streptococcus group A"),
                  coli = "-",       # Colistin
                  stringsAsFactors = FALSE)
  b <- data.frame(mo = c("Staphylococcus aureus",
                         "Streptococcus group A"),
                  coli = "R",       # Colistin
                  stringsAsFactors = FALSE)
  expect_equal(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)

  # piperacillin must be R in Enterobacteriaceae when tica is R
  library(dplyr)
  expect_equal(suppressWarnings(
    septic_patients %>%
      mutate(tica = as.rsi("R"),
             pipe = as.rsi("S")) %>%
      eucast_rules(col_mo = "mo") %>%
      left_join_microorganisms() %>%
      filter(family == "Enterobacteriaceae") %>%
      pull(pipe) %>%
      unique() %>%
      as.character()),
    "R")

  # azit and clar must be equal to eryt
  a <- suppressWarnings(
    septic_patients %>%
      transmute(mo,
                eryt,
                azit = as.rsi("R"),
                clar = as.rsi("R")) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(clar))
  b <-   suppressWarnings(
    septic_patients %>%
      select(mo, eryt) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(eryt))

  expect_identical(a[!is.na(b)],
                   b[!is.na(b)])

  # amox is inferred by benzylpenicillin in Kingella kingae
  expect_equal(
    suppressWarnings(
      as.list(eucast_rules(
        data.frame(mo = as.mo("Kingella kingae"),
                   peni = "S",
                   amox = "-",
                   stringsAsFactors = FALSE)
        , info = FALSE))$amox
    ),
    "S")

  # also test norf
  expect_output(suppressWarnings(eucast_rules(septic_patients %>% mutate(norf = "S", nali = "S"))))

  # check verbose output
  expect_output(suppressWarnings(eucast_rules(septic_patients, verbose = TRUE)))

})

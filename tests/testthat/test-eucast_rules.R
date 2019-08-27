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

context("eucast_rules.R")

test_that("EUCAST rules work", {

  # thoroughly check input table
  expect_equal(colnames(eucast_rules_file),
               c("if_mo_property", "like.is.one_of", "this_value",
                 "and_these_antibiotics", "have_these_values",
                 "then_change_these_antibiotics", "to_value",
                 "reference.rule", "reference.rule_group"))

  expect_error(suppressWarnings(eucast_rules(example_isolates, col_mo = "Non-existing")))
  expect_error(eucast_rules(x = "text"))
  expect_error(eucast_rules(data.frame(a = "test")))
  expect_error(eucast_rules(data.frame(mo = "test"), rules = "invalid rules set"))
  
  expect_warning(eucast_rules(data.frame(mo = "Escherichia coli", vancomycin = "S")))

  expect_identical(colnames(example_isolates),
                   colnames(suppressWarnings(eucast_rules(example_isolates))))

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

  a <- data.frame(mo = c("Staphylococcus aureus",
                         "Streptococcus group A"),
                  COL = "-",       # Colistin
                  stringsAsFactors = FALSE)
  b <- data.frame(mo = c("Staphylococcus aureus",
                         "Streptococcus group A"),
                  COL = "R",       # Colistin
                  stringsAsFactors = FALSE)
  expect_equal(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)

  # piperacillin must be R in Enterobacteriaceae when tica is R
  library(dplyr)
  expect_equal(suppressWarnings(
    example_isolates %>%
      mutate(TIC = as.rsi("R"),
             PIP = as.rsi("S")) %>%
      eucast_rules(col_mo = "mo") %>%
      left_join_microorganisms() %>%
      filter(family == "Enterobacteriaceae") %>%
      pull(PIP) %>%
      unique() %>%
      as.character()),
    "R")

  # Azithromicin and Clarythromycin must be equal to Erythromycin
  a <- suppressWarnings(
    example_isolates %>%
      transmute(mo,
                ERY,
                AZM = as.rsi("R"),
                CLR = as.rsi("R")) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(CLR))
  b <-   suppressWarnings(
    example_isolates %>%
      select(mo, ERY) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(ERY))

  expect_identical(a[!is.na(b)],
                   b[!is.na(b)])

  # amox is inferred by benzylpenicillin in Kingella kingae
  expect_equal(
    suppressWarnings(
      as.list(eucast_rules(
        data.frame(mo = as.mo("Kingella kingae"),
                   PEN = "S",
                   AMX = "-",
                   stringsAsFactors = FALSE)
        , info = FALSE))$AMX
    ),
    "S")

  # also test norf
  expect_output(suppressWarnings(eucast_rules(example_isolates %>% mutate(NOR = "S", NAL = "S"))))

  # check verbose output
  expect_output(suppressWarnings(eucast_rules(example_isolates, verbose = TRUE)))

})

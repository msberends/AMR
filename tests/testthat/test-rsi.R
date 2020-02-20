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

context("rsi.R")

test_that("rsi works", {
  
  skip_on_cran()
  
  expect_true(as.rsi("S") < as.rsi("I"))
  expect_true(as.rsi("I") < as.rsi("R"))
  expect_true(as.rsi("R") > as.rsi("S"))
  expect_true(is.rsi(as.rsi("S")))

  # print plots, should not raise errors
  barplot(as.rsi(c("S", "I", "R")))
  plot(as.rsi(c("S", "I", "R")))
  print(as.rsi(c("S", "I", "R")))

  expect_equal(suppressWarnings(as.logical(as.rsi("INVALID VALUE"))), NA)

  expect_equal(summary(as.rsi(c("S", "R"))), c("Class" = "rsi",
                                               "<NA>" = "0",
                                               "Sum S" = "1",
                                               "Sum IR" = "1",
                                               "-Sum R" = "1",
                                               "-Sum I" = "0"))

  expect_identical(as.logical(lapply(example_isolates, is.rsi.eligible)),
                   rep(FALSE, length(example_isolates)))

  library(dplyr)
  # 40 rsi columns
  expect_equal(example_isolates %>%
                 mutate_at(vars(PEN:RIF), as.character) %>%
                 lapply(is.rsi.eligible) %>%
                 as.logical() %>%
                 sum(),
               40)

})

test_that("mic2rsi works", {
  
  skip_on_cran()
  
  expect_equal(as.character(
    as.rsi(x = as.mic(0.125),
                      mo = "B_STRPT_PNMN",
                      ab = "AMX",
                      guideline = "EUCAST")),
    "S")
  expect_equal(as.character(
    as.rsi(x = as.mic(4),
           mo = "B_STRPT_PNMN",
           ab = "AMX",
           guideline = "EUCAST")),
    "R")

  expect_true(example_isolates %>%
                mutate(amox_mic = as.mic(2)) %>%
                select(mo, amox_mic) %>%
                as.rsi() %>%
                pull(amox_mic) %>%
                is.rsi())
  
  expect_warning(data.frame(mo = "E. coli",
                            NIT = c("<= 2", 32)) %>%
                   as.rsi())
  expect_message(data.frame(mo = "E. coli",
                            NIT = c("<= 2", 32),
                            uti = TRUE) %>%
                   as.rsi())
  expect_message(
    data.frame(mo = "E. coli",
               NIT = c("<= 2", 32),
               specimen = c("urine", "blood")) %>%
      as.rsi())
})

test_that("disk2rsi works", {
  
  skip_on_cran()
  
  expect_equal(as.character(
    as.rsi(x = as.disk(22),
           mo = "B_STRPT_PNMN",
           ab = "ERY",
           guideline = "CLSI")),
    "S")
  expect_equal(as.character(
    as.rsi(x = as.disk(18),
           mo = "B_STRPT_PNMN",
           ab = "ERY",
           guideline = "CLSI")),
    "I")
  expect_equal(as.character(
    as.rsi(x = as.disk(10),
           mo = "B_STRPT_PNMN",
           ab = "ERY",
           guideline = "CLSI")),
    "R")

  expect_true(example_isolates %>%
                mutate(amox_disk = as.disk(15)) %>%
                select(mo, amox_disk) %>%
                as.rsi(guideline = "CLSI") %>%
                pull(amox_disk) %>%
                is.rsi())
})

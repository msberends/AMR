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

context("filter_ab_class.R")

test_that("ATC-group filtering works", {
  library(dplyr)
  expect_gt(septic_patients %>% filter_ab_class("carbapenem") %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_aminoglycosides() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_carbapenems() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_cephalosporins() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_1st_cephalosporins() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_2nd_cephalosporins() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_3rd_cephalosporins() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_4th_cephalosporins() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_fluoroquinolones() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_glycopeptides() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_macrolides() %>% nrow(), 0)
  expect_gt(septic_patients %>% filter_tetracyclines() %>% nrow(), 0)

  expect_gt(septic_patients %>% filter_carbapenems("R", "all") %>% nrow(), 0)

  expect_error(septic_patients %>% filter_carbapenems(result = "test"))
  expect_error(septic_patients %>% filter_carbapenems(scope = "test"))
  expect_warning(septic_patients %>% select(1:3) %>% filter_carbapenems())
})

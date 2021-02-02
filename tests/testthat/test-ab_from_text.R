# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

context("ab_from_text.R")

test_that("ab_from_text works", {
  skip_on_cran()
  
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds")[[1]],
                   as.ab("Amoxicillin"))
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds", thorough_search = TRUE)[[1]],
                   as.ab("Amoxicillin"))
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds", thorough_search = FALSE)[[1]],
                   as.ab("Amoxicillin"))
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds", translate_ab = TRUE)[[1]],
                   "Amoxicillin")
  expect_identical(ab_from_text("administered amoxi/clav and cipro", collapse = ", ")[[1]],
                   "AMC, CIP")
  
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds", type = "dose")[[1]],
                   500)
  expect_identical(ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds", type = "admin")[[1]],
                   "oral")
})

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

context("pca.R")

test_that("PCA works", {
  library(dplyr)
  resistance_data <- example_isolates %>% 
    filter(mo %in% as.mo(c("E. coli", "K. pneumoniae", "S. aureus"))) %>% 
    select(mo, AMC, CXM, CTX, TOB, TMP) %>% 
    group_by(order = mo_order(mo),       # group on anything, like order
             genus = mo_genus(mo)) %>%   #  and genus as we do here
    summarise_if(is.rsi, resistance, minimum = 0)
  
  expect_s3_class(pca(resistance_data), "prcomp")
  expect_s3_class(prcomp(resistance_data), "prcomp")
  
  ggplot_pca(pca(resistance_data), ellipse = TRUE)
})

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
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
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
# ==================================================================== #

context("zzz.R")

test_that("imports work", {
  skip_on_cran()
  
  # Check if these function still exist in the package (all are in Suggests field)
  # Since GitHub Action runs every night, we will be emailed when a dependency fails based on this unit test
  import_functions <- c(
    "anti_join" = "dplyr",
    "cur_column" = "dplyr",
    "cur_data" = "dplyr",
    "document_position" = "rstudioapi",
    "document_range" = "rstudioapi",
    "freq.default" = "cleaner",
    "full_join" = "dplyr",
    "getSourceEditorContext" = "rstudioapi",
    "has_internet" = "curl",
    "html_attr" = "rvest",
    "html_children" = "rvest",
    "html_node" = "rvest",
    "html_nodes" = "rvest",
    "html_table" = "rvest",
    "html_text" = "rvest",
    "inline_hist" = "skimr",
    "inner_join" = "dplyr",
    "insertText" = "rstudioapi",
    "left_join" = "dplyr",
    "new_pillar_shaft_simple" = "pillar",
    "peek_vars" = "tidyselect",
    "read_excel" = "readxl",
    "read_html" = "xml2",
    "right_join" = "dplyr",
    "semi_join" = "dplyr",
    "sfl" = "skimr",
    "showQuestion" = "rstudioapi")
  
  for (i in seq_len(length(import_functions))) {
    fn <- names(import_functions)[i]
    pkg <- unname(import_functions[i])
    expect(!is.null(import_fn(name = fn, pkg = pkg, error_on_fail = FALSE)),
           failure_message = paste0("Function ", pkg, "::", fn, "() does not exist"))
  }
})

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

test_that("imports work", {
  skip_on_cran()
  
  import_functions <- c(
    cleaner = "freq.default",
    curl = "has_internet",
    dplyr = "cur_column",
    dplyr = "peek_mask",
    readxl = "read_excel",
    rstudioapi = "showQuestion",
    rvest = "html_attr",
    rvest = "html_children",
    rvest = "html_node",
    rvest = "html_nodes",
    rvest = "html_table",
    rvest = "html_text",
    tidyselect = "peek_vars",
    xml2 = "read_html")
  
  for (i in seq_len(length(import_functions))) {
    pkg <- names(import_functions)[i]
    fn <- unname(import_functions[i])
    expect(!is.null(import_fn(name = fn, pkg = pkg, error_on_fail = FALSE)),
           failure_message = paste0("Function ", pkg, "::", fn, "() does not exist"))
  }
})

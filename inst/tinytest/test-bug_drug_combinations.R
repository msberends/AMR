# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

b <- suppressWarnings(bug_drug_combinations(example_isolates))
expect_inherits(b, "bug_drug_combinations")
expect_stdout(suppressMessages(print(b)))
expect_true(is.data.frame(format(b)))
expect_true(is.data.frame(format(b, combine_IR = TRUE, add_ab_group = FALSE)))
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  expect_true(example_isolates %>% 
                group_by(hospital_id) %>% 
                bug_drug_combinations(FUN = mo_gramstain) %>% 
                is.data.frame())
}

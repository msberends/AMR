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

# context("All examples")
# 
# # run all examples (will take forever)
# exported_functions <- ls("package:AMR")
# 
# for (i in seq_len(length(exported_functions))) {
#   test_that(paste(exported_functions[i], "works"), {
#     skip_on_cran()
#     expect_output(suppressWarnings(example(exported_functions[i], 
#                                            package = "AMR", 
#                                            give.lines = TRUE,
#                                            run.dontrun = TRUE,
#                                            run.donttest = TRUE)),
#                   label = paste0("Examples of function ", exported_functions[i]))
#   })
# }

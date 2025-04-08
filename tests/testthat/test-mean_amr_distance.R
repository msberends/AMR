# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

test_that("test-mean_amr_distance.R", {
  skip_on_cran()

  vctr_disk <- as.disk(c(20:25))
  vctr_mic <- as.mic(2^c(0:5))
  vctr_sir <- as.sir(c("S", "S", "I", "I", "R", "R"))

  expect_identical(
    mean_amr_distance(vctr_disk),
    (as.double(vctr_disk) - mean(as.double(vctr_disk))) / sd(as.double(vctr_disk))
  )

  expect_identical(
    mean_amr_distance(vctr_mic),
    (log2(vctr_mic) - mean(log2(vctr_mic))) / sd(log2(vctr_mic))
  )

  expect_identical(
    mean_amr_distance(vctr_sir, combine_SI = FALSE),
    (c(1, 1, 2, 2, 3, 3) - mean(c(1, 1, 2, 2, 3, 3))) / sd(c(1, 1, 2, 2, 3, 3))
  )
  expect_identical(
    mean_amr_distance(vctr_sir, combine_SI = TRUE),
    (c(1, 1, 1, 1, 3, 3) - mean(c(1, 1, 1, 1, 3, 3))) / sd(c(1, 1, 1, 1, 3, 3))
  )

  expect_equal(
    mean_amr_distance(data.frame(AMX = vctr_mic, GEN = vctr_sir, TOB = vctr_disk)),
    c(-1.10603655, -0.74968823, -0.39333990, -0.03699158, 0.96485397, 1.32120229),
    tolerance = 0.00001
  )

  expect_equal(
    mean_amr_distance(data.frame(AMX = vctr_mic, GEN = vctr_sir, TOB = vctr_disk), 2:3),
    c(-0.9909017, -0.7236405, -0.4563792, -0.1891180, 1.0463891, 1.3136503),
    tolerance = 0.00001
  )
})

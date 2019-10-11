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

context("get_locale.R")

test_that("get_locale works", {
  expect_identical(mo_genus("B_GRAMP", language = "pt"),
                   "(Gram positivos desconhecidos)")

  expect_identical(mo_fullname("CoNS", "en"), "Coagulase-negative Staphylococcus (CoNS)")
  expect_identical(mo_fullname("CoNS", "de"), "Koagulase-negative Staphylococcus (KNS)")
  expect_identical(mo_fullname("CoNS", "nl"), "Coagulase-negatieve Staphylococcus (CNS)")
  expect_identical(mo_fullname("CoNS", "es"), "Staphylococcus coagulasa negativo (SCN)")
  expect_identical(mo_fullname("CoNS", "it"), "Staphylococcus negativo coagulasi (CoNS)")
  expect_identical(mo_fullname("CoNS", "pt"), "Staphylococcus coagulase negativo (CoNS)")

})

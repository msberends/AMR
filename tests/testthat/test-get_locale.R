# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

context("get_locale.R")

test_that("get_locale works", {
  expect_identical(mo_genus("B_GRAMP", language = "pt"),
                   "(Gram positivos desconhecidos)")

  expect_identical(mo_fullname("CoNS", "en"), "Coagulase Negative Staphylococcus (CoNS)")
  expect_identical(mo_fullname("CoNS", "de"), "Koagulase-negative Staphylococcus (KNS)")
  expect_identical(mo_fullname("CoNS", "nl"), "Coagulase-negatieve Staphylococcus (CNS)")
  expect_identical(mo_fullname("CoNS", "es"), "Staphylococcus coagulasa negativo (CoNS)")
  expect_identical(mo_fullname("CoNS", "it"), "Staphylococcus negativo coagulasi (CoNS)")
  # expect_identical(mo_fullname("CoNS", "fr"), "Staphylococcus \u00e0 coagulase n\u00e9gative (CoNS)")
  expect_identical(mo_fullname("CoNS", "pt"), "Staphylococcus coagulase negativo (CoNS)")

})

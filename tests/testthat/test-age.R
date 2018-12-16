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

context("age.R")

test_that("age works", {
  expect_equal(age(x = c("1980-01-01", "1985-01-01", "1990-01-01"),
                   y = "2019-01-01"),
               c(39, 34, 29))

  expect_error(age(x = c("1980-01-01", "1985-01-01", "1990-01-01"),
                   y = c("2019-01-01", "2019-01-01")))

  expect_error(age(x = c("1980-01-01", "1985-01-01", "1990-01-01"),
                   y = "1975-01-01"))

  expect_warning(age(x = c("1800-01-01", "1805-01-01", "1810-01-01"),
                     y = "2019-01-01"))
})

test_that("age_groups works", {
  ages <- c(3, 8, 16, 54, 31, 76, 101, 43, 21)

  expect_equal(length(unique(age_groups(ages, 50))),
               2)
  expect_equal(length(unique(age_groups(ages, c(50, 60)))),
               3)
  expect_identical(class(age_groups(ages, "child")),
                   c("ordered", "factor"))

  expect_identical(class(age_groups(ages, "elderly")),
                   c("ordered", "factor"))

  expect_identical(class(age_groups(ages, "tens")),
                   c("ordered", "factor"))

  expect_identical(class(age_groups(ages, "fives")),
                   c("ordered", "factor"))

})

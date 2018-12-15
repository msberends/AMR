context("g.test.R")

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
  expect_identical(class(age_groups(ages)),
                   c("ordered", "factor"))


})

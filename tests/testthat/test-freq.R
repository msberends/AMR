context("freq.R")

test_that("frequency table works", {
  expect_equal(nrow(freq(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))), 5)
  expect_equal(nrow(frequency_tbl(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))), 5)

  # date column of septic_patients should contain 1151 unique dates
  expect_equal(nrow(freq(septic_patients$date)), 1151)
  expect_equal(nrow(freq(septic_patients$date)),
               length(unique(septic_patients$date)))

  # character
  expect_output(print(freq(septic_patients$bactid)))
  # integer
  expect_output(print(freq(septic_patients$age)))
  # date
  expect_output(print(freq(septic_patients$date)))
  # factor
  expect_output(print(freq(septic_patients$hospital_id)))
  # table
  expect_output(print(freq(table(septic_patients$sex, septic_patients$age))))

  library(dplyr)
  expect_output(septic_patients %>% select(1:2) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:3) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:4) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:5) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:6) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:7) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:8) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:9) %>% freq() %>% print())

  # top 5
  expect_equal(
    septic_patients %>%
      freq(bactid) %>%
      top_freq(5) %>%
      length(),
    5)
  # there're more than 5 lowest values
  expect_gt(
    septic_patients %>%
      freq(bactid) %>%
      top_freq(-5) %>%
      length(),
    5)
  # n has length > 1
  expect_error(
    septic_patients %>%
      freq(bactid) %>%
      top_freq(n = c(1, 2))
  )
  # input must be freq tbl
  expect_error(septic_patients %>% top_freq(1))

  # charts from plot and hist, should not raise errors
  plot(freq(septic_patients, age))
  hist(freq(septic_patients, age))

  # check vector
  expect_identical(septic_patients %>%
                     freq(age) %>%
                     as.vector() %>%
                     sort(),
                   septic_patients %>%
                     pull(age) %>%
                     sort())

  # check format
  expect_identical(septic_patients %>%
                     freq(age) %>%
                     format() %>%
                     apply(2, class) %>%
                     unname(),
                   rep("character", 5))

  # check tibble
  expect_identical(septic_patients %>%
                     freq(age) %>%
                     as_tibble() %>%
                     class() %>%
                     .[1],
                   "tbl_df")
})


context("freq.R")

test_that("frequency table works", {
  library(dplyr)

  expect_equal(nrow(freq(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))), 5)
  expect_equal(nrow(frequency_tbl(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))), 5)

  # date column of septic_patients should contain 1151 unique dates
  expect_equal(nrow(freq(septic_patients$date)), 1151)
  expect_equal(nrow(freq(septic_patients$date)),
               length(unique(septic_patients$date)))

  expect_output(print(septic_patients %>% freq(age, nmax = Inf)))
  expect_output(print(freq(septic_patients$age, nmax = Inf)))
  expect_output(print(freq(septic_patients$age, nmax = NA)))
  expect_output(print(freq(septic_patients$age, nmax = NULL)))
  expect_output(print(freq(septic_patients$age, sort.count = FALSE)))
  expect_output(print(freq(septic_patients$age, markdown = TRUE)))
  expect_output(print(freq(septic_patients$age, markdown = TRUE), markdown = FALSE))
  expect_output(print(freq(septic_patients$age, markdown = TRUE), markdown = TRUE))
  expect_output(print(freq(septic_patients$age[0])))

  # character
  expect_output(print(freq(septic_patients$mo)))
  # integer
  expect_output(print(freq(septic_patients$age)))
  # date
  expect_output(print(freq(septic_patients$date)))
  # factor
  expect_output(print(freq(septic_patients$hospital_id)))
  # table
  expect_output(print(freq(table(septic_patients$sex, septic_patients$age))))
  # rsi
  expect_output(print(freq(septic_patients$amcl)))
  # hms
  expect_output(suppressWarnings(print(freq(hms::as.hms(sample(c(0:86399), 50))))))
  # matrix
  expect_output(print(freq(as.matrix(septic_patients$age))))
  expect_output(print(freq(as.matrix(septic_patients[, c("age", "sex")]))))
  # list
  expect_output(print(freq(list(age = septic_patients$age))))
  expect_output(print(freq(list(age = septic_patients$age, sex = septic_patients$sex))))

  library(dplyr)
  expect_output(septic_patients %>% select(1:2) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:3) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:4) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:5) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:6) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:7) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:8) %>% freq() %>% print())
  expect_output(septic_patients %>% select(1:9) %>% freq() %>% print())
  expect_output(print(freq(septic_patients$age), nmax = 20))

  # top 5
  expect_equal(
    septic_patients %>%
      freq(mo) %>%
      top_freq(5) %>%
      length(),
    5)
  # there're more than 5 lowest values
  expect_gt(
    septic_patients %>%
      freq(mo) %>%
      top_freq(-5) %>%
      length(),
    5)
  # n has length > 1
  expect_error(
    septic_patients %>%
      freq(mo) %>%
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

  expect_error(septic_patients %>% freq(nonexisting))
  expect_error(septic_patients %>% select(1:10) %>% freq())
  expect_error(septic_patients %>% freq(peni, oxac, clox, amox, amcl,
                                        ampi, pita, czol, cfep, cfur))

})


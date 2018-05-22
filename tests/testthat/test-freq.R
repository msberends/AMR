context("freq.R")

test_that("frequency table works", {
  expect_equal(nrow(freq(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5), as.data.frame = TRUE)), 5)
  expect_equal(nrow(frequency_tbl(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5), as.data.frame = TRUE)), 5)

  # date column of septic_patients should contain 1662 unique dates
  expect_equal(nrow(freq(septic_patients$date, as.data.frame = TRUE)), 1662)
  expect_equal(nrow(freq(septic_patients$date, as.data.frame = TRUE)),
               length(unique(septic_patients$date)))

  expect_output(freq(septic_patients$age))
  expect_output(freq(septic_patients$date))
  expect_output(freq(septic_patients$hospital_id))

  library(dplyr)
  expect_output(septic_patients %>% select(1:2) %>% freq())
  expect_output(septic_patients %>% select(1:3) %>% freq())
  expect_output(septic_patients %>% select(1:4) %>% freq())
  expect_output(septic_patients %>% select(1:5) %>% freq())
  expect_output(septic_patients %>% select(1:6) %>% freq())
  expect_output(septic_patients %>% select(1:7) %>% freq())
  expect_output(septic_patients %>% select(1:8) %>% freq())
  expect_output(septic_patients %>% select(1:9) %>% freq())
})


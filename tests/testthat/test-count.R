context("count.R")

test_that("counts work", {
  # amox resistance in `septic_patients`
  expect_equal(count_R(septic_patients$amox), 662)
  expect_equal(count_I(septic_patients$amox), 3)
  expect_equal(count_S(septic_patients$amox), 335)
  expect_equal(count_R(septic_patients$amox) + count_I(septic_patients$amox),
               count_IR(septic_patients$amox))
  expect_equal(count_S(septic_patients$amox) + count_I(septic_patients$amox),
               count_SI(septic_patients$amox))

  library(dplyr)
  expect_equal(septic_patients %>% count_S(amcl), 1057)
  expect_equal(septic_patients %>% count_S(amcl, gent), 1396)
  expect_equal(septic_patients %>% count_all(amcl, gent), 1517)
  expect_identical(septic_patients %>% count_all(amcl, gent),
                   septic_patients %>% count_S(amcl, gent) +
                     septic_patients %>% count_IR(amcl, gent))

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro = count_S(cipr),
                           genta = count_S(gent),
                           combination = count_S(cipr, gent)) %>%
                 pull(combination),
               c(192, 446, 184, 474))

  # count_df
  expect_equal(
    septic_patients %>% select(amox) %>% count_df() %>% pull(Value),
    c(septic_patients$amox %>% count_S(),
      septic_patients$amox %>% count_I(),
      septic_patients$amox %>% count_R())
  )
  expect_equal(
    septic_patients %>% select(amox) %>% count_df(combine_IR = TRUE) %>% pull(Value),
    c(septic_patients$amox %>% count_S(),
      septic_patients$amox %>% count_IR())
  )

  # warning for speed loss
  expect_warning(count_R(as.character(septic_patients$amcl)))
  expect_warning(count_I(as.character(septic_patients$amcl)))
  expect_warning(count_S(as.character(septic_patients$amcl,
                                      septic_patients$gent)))
  expect_warning(count_S(septic_patients$amcl,
                         as.character(septic_patients$gent)))

  # check for errors
  expect_error(count_IR("test", minimum = "test"))
  expect_error(count_IR("test", as_percent = "test"))
  expect_error(count_I("test", minimum = "test"))
  expect_error(count_I("test", as_percent = "test"))
  expect_error(count_S("test", minimum = "test"))
  expect_error(count_S("test", as_percent = "test"))

  expect_error(count_df(c("A", "B", "C")))
  expect_error(count_df(septic_patients[,"date"]))

})

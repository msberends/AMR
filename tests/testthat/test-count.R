context("count.R")

test_that("counts work", {
  # amox resistance in `septic_patients`
  expect_equal(count_R(septic_patients$amox), 659)
  expect_equal(count_I(septic_patients$amox), 3)
  expect_equal(count_S(septic_patients$amox), 336)
  expect_equal(count_R(septic_patients$amox) + count_I(septic_patients$amox),
               count_IR(septic_patients$amox))
  expect_equal(count_S(septic_patients$amox) + count_I(septic_patients$amox),
               count_SI(septic_patients$amox))

  expect_equal(septic_patients %>% count_S(amcl), 1056)
  expect_equal(septic_patients %>% count_S(amcl, gent), 1385)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro = count_S(cipr),
                           genta = count_S(gent),
                           combination = count_S(cipr, gent)) %>%
                 pull(combination),
               c(192, 440, 184, 474))

  expect_equal(septic_patients %>% select(amox, cipr) %>% count_df(translate_ab = "official") %>% nrow(),
               6)

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

})

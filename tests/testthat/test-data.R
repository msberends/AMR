context("data.R")

test_that("data sets are valid", {
  # IDs should always be unique
  expect_identical(nrow(antibiotics), length(unique(antibiotics$atc)))
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))

  # there should be no diacritics (i.e. non ASCII) characters in the datasets
  library(dplyr)
  # check only character variables:
  test_microorganisms <- microorganisms %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_microorganisms.old <- microorganisms.old %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_antibiotics <- antibiotics %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  test_septic_patients <- septic_patients %>% select_if(is.character) %>% as.data.frame(stringsAsFactors = FALSE)
  # and compare them with their transformed version:
  expect_identical(test_microorganisms,
                   test_microorganisms %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_microorganisms.old,
                   test_microorganisms.old %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_antibiotics,
                   test_antibiotics %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

  expect_identical(test_septic_patients,
                   test_septic_patients %>%
                     lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
                     as.data.frame(stringsAsFactors = FALSE))

})

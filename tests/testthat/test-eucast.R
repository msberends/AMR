context("eucast.R")

test_that("EUCAST rules work", {

  expect_error(EUCAST_rules(septic_patients, col_bactid = "Non-existing"))


  expect_identical(colnames(septic_patients),
                   colnames(suppressWarnings(EUCAST_rules(septic_patients))))

  a <- data.frame(bactid =
                    c("KLEPNE",  # Klebsiella pneumoniae
                      "PSEAER",  # Pseudomonas aeruginosa
                      "ENTAER"), # Enterobacter aerogenes
                  amox = "-",           # Amoxicillin
                  stringsAsFactors = FALSE)
  b <- data.frame(bactid =
                    c("KLEPNE",  # Klebsiella pneumoniae
                      "PSEAER",  # Pseudomonas aeruginosa
                      "ENTAER"), # Enterobacter aerogenes
                  amox = "R",           # Amoxicillin
                  stringsAsFactors = FALSE)
  expect_warning(EUCAST_rules(a, info = FALSE))
  expect_identical(suppressWarnings(EUCAST_rules(a, info = FALSE)), b)
  expect_identical(suppressWarnings(interpretive_reading(a, info = TRUE)), b)

  a <- data.frame(bactid =
                    c("STAAUR",  # Staphylococcus aureus
                      "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "-",           # Colistin
                  stringsAsFactors = FALSE)
  b <- data.frame(bactid =
                    c("STAAUR",  # Staphylococcus aureus
                      "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "R",           # Colistin
                  stringsAsFactors = FALSE)
  expect_equal(suppressWarnings(EUCAST_rules(a, info = FALSE)), b)

  # pita must be R in Enterobacteriaceae when tica is R
  library(dplyr)
  expect_equal(suppressWarnings(
    septic_patients %>%
      mutate(tica = as.rsi("R"),
             pita = as.rsi("S")) %>%
      EUCAST_rules(col_bactid = "bactid") %>%
      left_join_microorganisms() %>%
      filter(family == "Enterobacteriaceae") %>%
      pull(pita) %>%
      unique() %>%
      as.character()),
    "R")
  # azit and clar must be equal to eryt
  expect_equal(suppressWarnings(
    septic_patients %>%
      mutate(azit = as.rsi("R"),
             clar = as.rsi("R")) %>%
      EUCAST_rules(col_bactid = "bactid") %>%
      pull(clar)),
    suppressWarnings(
      septic_patients %>%
        EUCAST_rules(col_bactid = "bactid") %>%
        left_join_microorganisms() %>%
        pull(eryt)))

})

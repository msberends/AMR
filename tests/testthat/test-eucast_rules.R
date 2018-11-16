context("eucast_rules.R")

test_that("EUCAST rules work", {

  expect_error(suppressWarnings(eucast_rules(septic_patients, col_mo = "Non-existing")))

  expect_identical(colnames(septic_patients),
                   colnames(suppressWarnings(eucast_rules(septic_patients))))

  a <- data.frame(mo = c("KLEPNE",  # Klebsiella pneumoniae
                         "PSEAER",  # Pseudomonas aeruginosa
                         "ENTAER"), # Enterobacter aerogenes
                  amox = "-",           # Amoxicillin
                  stringsAsFactors = FALSE)
  b <- data.frame(mo = c("KLEPNE",  # Klebsiella pneumoniae
                         "PSEAER",  # Pseudomonas aeruginosa
                         "ENTAER"), # Enterobacter aerogenes
                  amox = "R",       # Amoxicillin
                  stringsAsFactors = FALSE)
  expect_identical(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)
  expect_identical(suppressWarnings(eucast_rules(a, "mo", info = TRUE)), b)
  expect_identical(suppressWarnings(interpretive_reading(a, "mo", info = TRUE)), b)

  a <- data.frame(mo = c("STAAUR",  # Staphylococcus aureus
                         "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "-",       # Colistin
                  stringsAsFactors = FALSE)
  b <- data.frame(mo = c("STAAUR",  # Staphylococcus aureus
                         "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "R",       # Colistin
                  stringsAsFactors = FALSE)
  expect_equal(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)

  # piperacillin must be R in Enterobacteriaceae when tica is R
  library(dplyr)
  expect_equal(suppressWarnings(
    septic_patients %>%
      mutate(tica = as.rsi("R"),
             pipe = as.rsi("S")) %>%
      eucast_rules(col_mo = "mo") %>%
      left_join_microorganisms() %>%
      filter(family == "Enterobacteriaceae") %>%
      pull(pipe) %>%
      unique() %>%
      as.character()),
    "R")

  # azit and clar must be equal to eryt
  a <- suppressWarnings(
    septic_patients %>%
      transmute(mo,
                eryt,
                azit = as.rsi("R"),
                clar = as.rsi("R")) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(clar))
  b <-   suppressWarnings(
    septic_patients %>%
      select(mo, eryt) %>%
      eucast_rules(col_mo = "mo") %>%
      pull(eryt))

  expect_identical(a[!is.na(b)],
                   b[!is.na(b)])

  # amox is inferred by benzylpenicillin in Kingella kingae
  expect_equal(
    as.list(eucast_rules(
      data.frame(mo = as.mo("Kingella kingae"),
                 peni = "S",
                 amox = "-",
                 stringsAsFactors = FALSE)
      , info = FALSE))$amox,
    "S")

  expect_output(suppressWarnings(eucast_rules(septic_patients, verbose = TRUE)))

})

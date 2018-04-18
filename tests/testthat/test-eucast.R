context("eucast.R")

test_that("EUCAST rules work", {
  a <- data.frame(bactid = c("KLEPNE",  # Klebsiella pneumoniae
                             "PSEAER",  # Pseudomonas aeruginosa
                             "ENTAER"), # Enterobacter aerogenes
                  amox = "-",           # Amoxicillin
                  stringsAsFactors = FALSE)
  b <- data.frame(bactid = c("KLEPNE",  # Klebsiella pneumoniae
                             "PSEAER",  # Pseudomonas aeruginosa
                             "ENTAER"), # Enterobacter aerogenes
                  amox = "R",           # Amoxicillin
                  stringsAsFactors = FALSE)
  expect_equal(EUCAST_rules(a, info = FALSE), b)
  expect_equal(suppressWarnings(interpretive_reading(a, info = TRUE)), b)

  a <- data.frame(bactid = c("STAAUR",  # Staphylococcus aureus
                             "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "-",           # Colistin
                  stringsAsFactors = FALSE)
  b <- data.frame(bactid = c("STAAUR",  # Staphylococcus aureus
                             "STCGRA"), # Streptococcus pyognenes (Lancefield Group A)
                  coli = "R",           # Colistin
                  stringsAsFactors = FALSE)
  expect_equal(EUCAST_rules(a, info = FALSE), b)
})

test_that("MO properties work", {
  expect_equal(mo_property("ESCCOL"), "Escherichia coli")
  expect_equal(mo_property("STAAUR"), "Staphylococcus aureus")
})

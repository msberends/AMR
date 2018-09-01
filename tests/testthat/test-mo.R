context("mo.R")

test_that("as.mo works", {
  expect_identical(
    as.character(as.mo(c("E. coli", "H. influenzae"))),
    c("ESCCOL", "HAEINF"))

  expect_equal(as.character(as.mo("Escherichia coli")), "ESCCOL")
  expect_equal(as.character(as.mo("Escherichia  coli")), "ESCCOL")
  expect_equal(as.character(as.mo("Escherichia  species")), "ESC")
  expect_equal(as.character(as.mo(" ESCCOL ")), "ESCCOL")
  expect_equal(as.character(as.mo("klpn")), "KLEPNE")
  expect_equal(as.character(as.mo("Klebsiella")), "KLE")
  expect_equal(as.character(as.mo("K. pneu rhino")), "KLEPNERH") # K. pneumoniae subspp. rhinoscleromatis
  expect_equal(as.character(as.mo("coagulase negative")), "STACNS")

  expect_equal(as.character(as.mo("P. aer")), "PSEAER") # not Pasteurella aerogenes

  expect_equal(as.character(as.mo("Negative rods")), "GNR")
  expect_equal(as.character(as.mo("Gram negative rods")), "GNR")

  # GLIMS
  expect_equal(as.character(as.mo("shiboy")), "SHIBOY")

  expect_equal(as.character(as.mo("MRSE")), "STAEPI")
  expect_equal(as.character(as.mo("VRE")), "ENC")
  expect_equal(as.character(as.mo("MRPA")), "PSEAER")
  expect_equal(as.character(as.mo("PISP")), "STCPNE")
  expect_equal(as.character(as.mo("PRSP")), "STCPNE")
  expect_equal(as.character(as.mo("VISP")), "STCPNE")
  expect_equal(as.character(as.mo("VRSP")), "STCPNE")

  expect_identical(
    as.character(
      as.mo(c("stau",
                     "STAU",
                     "staaur",
                     "S. aureus",
                     "S aureus",
                     "Staphylococcus aureus",
                     "MRSA",
                     "VISA"))),
    rep("STAAUR", 8))

  # check for Becker classification
  expect_identical(as.character(guess_mo("S. epidermidis", Becker = FALSE)), "STAEPI")
  expect_identical(as.character(guess_mo("S. epidermidis", Becker = TRUE)),  "STACNS")
  expect_identical(as.character(guess_mo("STAEPI",         Becker = TRUE)),  "STACNS")
  expect_identical(as.character(guess_mo("S. intermedius", Becker = FALSE)), "STAINT")
  expect_identical(as.character(guess_mo("S. intermedius", Becker = TRUE)),  "STACPS")
  expect_identical(as.character(guess_mo("STAINT",         Becker = TRUE)),  "STACPS")
  # aureus must only be influenced if Becker = "all"
  expect_identical(as.character(guess_mo("STAAUR", Becker = FALSE)), "STAAUR")
  expect_identical(as.character(guess_mo("STAAUR", Becker = TRUE)),  "STAAUR")
  expect_identical(as.character(guess_mo("STAAUR", Becker = "all")), "STACPS")

  # check for Lancefield classification
  expect_identical(as.character(guess_mo("S. pyogenes", Lancefield = FALSE)), "STCPYO")
  expect_identical(as.character(guess_mo("S. pyogenes", Lancefield = TRUE)),  "STCGRA")
  expect_identical(as.character(guess_mo("STCPYO",      Lancefield = TRUE)),  "STCGRA")
  expect_identical(as.character(guess_mo("S. agalactiae",  Lancefield = FALSE)),  "STCAGA")
  expect_identical(as.character(guess_mo("S. agalactiae",  Lancefield = TRUE)),   "STCGRB") # group B
  expect_identical(as.character(guess_mo("S. equisimilis", Lancefield = FALSE)),  "STCEQS")
  expect_identical(as.character(guess_mo("S. equisimilis", Lancefield = TRUE)),   "STCGRC") # group C
  expect_identical(as.character(guess_mo("S. anginosus",   Lancefield = FALSE)),  "STCANG")
  expect_identical(as.character(guess_mo("S. anginosus",   Lancefield = TRUE)),   "STCGRF") # group F
  expect_identical(as.character(guess_mo("S. sanguis",     Lancefield = FALSE)),  "STCSAN")
  expect_identical(as.character(guess_mo("S. sanguis",     Lancefield = TRUE)),   "STCGRH") # group H
  expect_identical(as.character(guess_mo("S. salivarius",  Lancefield = FALSE)),  "STCSAL")
  expect_identical(as.character(guess_mo("S. salivarius",  Lancefield = TRUE)),   "STCGRK") # group K

  library(dplyr)

  # select with one column
  expect_identical(
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus) %>%
      as.mo() %>%
      as.character(),
    c("ESC", "ESC", "STA", "STA", "STA",
      "STA", "STA", "STA", "STA", "STA"))

  # select with two columns
  expect_identical(
    septic_patients[1:10,] %>%
      pull(mo),
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus, species) %>%
      as.mo() %>%
      as.character())

  # unknown results
  expect_warning(as.mo(c("INVALID", "Yeah, unknown")))

  # too many columns
  expect_error(septic_patients %>% select(1:3) %>% as.mo())

  # print
  expect_output(print(as.mo(c("ESCCOL", NA))))

  # helper function
  expect_identical(as.mo("ESCCOL"),
                   guess_mo("ESCCOL"))

  # test pull
  expect_equal(nrow(septic_patients %>% mutate(mo = as.mo(mo))),
               2000)

  # test data.frame
  expect_equal(nrow(data.frame(test = as.mo("ESCCOL"))),
               1)

  # check empty values
  expect_equal(as.character(suppressWarnings(as.mo(""))),
               NA_character_)

})

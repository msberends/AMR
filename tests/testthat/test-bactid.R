context("bactid.R")

test_that("as.bactid works", {
  expect_identical(
    as.character(as.bactid(c("E. coli", "H. influenzae"))),
    c("ESCCOL", "HAEINF"))

  expect_equal(as.character(as.bactid("Escherichia coli")), "ESCCOL")
  expect_equal(as.character(as.bactid("Escherichia  coli")), "ESCCOL")
  expect_equal(as.character(as.bactid("Escherichia  species")), "ESC")
  expect_equal(as.character(as.bactid(" ESCCOL ")), "ESCCOL")
  expect_equal(as.character(as.bactid("klpn")), "KLEPNE")

  expect_equal(as.character(as.bactid("P. aer")), "PSEAER") # not Pasteurella aerogenes

  expect_equal(as.character(as.bactid("Negative rods")), "GNR")

  # GLIMS
  expect_equal(as.character(as.bactid("shiboy")), "SHIBOY")

  expect_equal(as.character(as.bactid("MRSE")), "STAEPI")
  expect_equal(as.character(as.bactid("VRE")), "ENC")
  expect_equal(as.character(as.bactid("MRPA")), "PSEAER")
  expect_equal(as.character(as.bactid("PISP")), "STCPNE")
  expect_equal(as.character(as.bactid("PRSP")), "STCPNE")
  expect_equal(as.character(as.bactid("VISP")), "STCPNE")
  expect_equal(as.character(as.bactid("VRSP")), "STCPNE")

  expect_identical(
    as.character(
      as.bactid(c("stau",
                     "STAU",
                     "staaur",
                     "S. aureus",
                     "S aureus",
                     "Staphylococcus aureus",
                     "MRSA",
                     "VISA"))),
    rep("STAAUR", 8))

  # check for Becker classification
  expect_identical(as.character(guess_bactid("S. epidermidis", Becker = FALSE)), "STAEPI")
  expect_identical(as.character(guess_bactid("S. epidermidis", Becker = TRUE)),  "STACNS")
  expect_identical(as.character(guess_bactid("STAEPI",         Becker = TRUE)),  "STACNS")
  expect_identical(as.character(guess_bactid("S. intermedius", Becker = FALSE)), "STAINT")
  expect_identical(as.character(guess_bactid("S. intermedius", Becker = TRUE)),  "STACPS")
  expect_identical(as.character(guess_bactid("STAINT",         Becker = TRUE)),  "STACPS")
  # aureus must only be influenced if Becker = "all"
  expect_identical(as.character(guess_bactid("STAAUR", Becker = FALSE)), "STAAUR")
  expect_identical(as.character(guess_bactid("STAAUR", Becker = TRUE)),  "STAAUR")
  expect_identical(as.character(guess_bactid("STAAUR", Becker = "all")), "STACPS")

  # check for Lancefield classification
  expect_identical(as.character(guess_bactid("S. pyogenes", Lancefield = FALSE)), "STCPYO")
  expect_identical(as.character(guess_bactid("S. pyogenes", Lancefield = TRUE)),  "STCGRA")
  expect_identical(as.character(guess_bactid("STCPYO",      Lancefield = TRUE)),  "STCGRA")
  expect_identical(as.character(guess_bactid("S. agalactiae",  Lancefield = FALSE)),  "STCAGA")
  expect_identical(as.character(guess_bactid("S. agalactiae",  Lancefield = TRUE)),   "STCGRB") # group B
  expect_identical(as.character(guess_bactid("S. equisimilis", Lancefield = FALSE)),  "STCEQS")
  expect_identical(as.character(guess_bactid("S. equisimilis", Lancefield = TRUE)),   "STCGRC") # group C
  expect_identical(as.character(guess_bactid("S. anginosus",   Lancefield = FALSE)),  "STCANG")
  expect_identical(as.character(guess_bactid("S. anginosus",   Lancefield = TRUE)),   "STCGRF") # group F
  expect_identical(as.character(guess_bactid("S. sanguis",     Lancefield = FALSE)),  "STCSAN")
  expect_identical(as.character(guess_bactid("S. sanguis",     Lancefield = TRUE)),   "STCGRH") # group H
  expect_identical(as.character(guess_bactid("S. salivarius",  Lancefield = FALSE)),  "STCSAL")
  expect_identical(as.character(guess_bactid("S. salivarius",  Lancefield = TRUE)),   "STCGRK") # group K

  # select with one column
  expect_identical(
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus) %>%
      as.bactid() %>%
      as.character(),
    c("ESC", "ESC", "STA", "STA", "STA",
      "STA", "STA", "STA", "STA", "STA"))

  # select with two columns
  expect_identical(
    septic_patients[1:10,] %>%
      pull(bactid),
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus, species) %>%
      as.bactid() %>%
      as.character())

  # unknown results
  expect_warning(as.bactid(c("INVALID", "Yeah, unknown")))

  # print
  expect_output(print(as.bactid(c("ESCCOL", NA))))

  # helper function
  expect_identical(as.bactid("ESCCOL"),
                   guess_bactid("ESCCOL"))

  # test pull
  expect_equal(nrow(septic_patients %>% mutate(bactid = as.bactid(bactid))),
               2000)

  # test data.frame
  expect_equal(nrow(data.frame(test = as.bactid("ESCCOL"))),
               1)

  # check empty values
  expect_equal(as.character(suppressWarnings(as.bactid(""))),
               NA_character_)

})

test_that("bactid.property works", {
  expect_equal(bactid.family("E. coli"), "Enterobacteriaceae")
  expect_equal(bactid.genus("E. coli"), "Escherichia")
  expect_equal(bactid.species("E. coli"), "coli")
  expect_equal(bactid.subspecies("E. coli"), NA_character_)
  expect_equal(bactid.fullname("E. coli"), "Escherichia coli")
  expect_equal(bactid.type("E. coli"), "Bacteria")
  expect_equal(bactid.gramstain("E. coli"), "Negative rods")
  expect_equal(bactid.aerobic("E. coli"), TRUE)
  expect_equal(bactid.type_nl("E. coli"), "Bacterie")
  expect_equal(bactid.gramstain_nl("E. coli"), "Negatieve staven")
})

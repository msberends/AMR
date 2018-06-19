context("guess_bactid.R")

test_that("guess_bactid works", {
  expect_identical(
    guess_bactid(c("E. coli", "H. influenzae")),
    c("ESCCOL", "HAEINF"))

  expect_equal(guess_bactid("Escherichia coli"), "ESCCOL")
  expect_equal(guess_bactid("P. aer"), "PSEAER") # not Pasteurella aerogenes

  expect_equal(guess_bactid("Negative rods"), "GNR")

  expect_equal(guess_bactid("MRSE"), "STAEPI")
  expect_equal(guess_bactid("VRE"), "ENC")
  expect_equal(guess_bactid("MRPA"), "PSEAER")
  expect_equal(guess_bactid("PISP"), "STCPNE")
  expect_equal(guess_bactid("PRSP"), "STCPNE")
  expect_equal(guess_bactid("VISP"), "STCPNE")
  expect_equal(guess_bactid("VRSP"), "STCPNE")

  expect_identical(
    guess_bactid(c("stau",
                   "STAU",
                   "staaur",
                   "S. aureus",
                   "S aureus",
                   "Staphylococcus aureus",
                   "MRSA",
                   "VISA")),
    rep("STAAUR", 8))

  # select with one column
  expect_identical(
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus) %>%
      guess_bactid(),
    c("STC", "STC", "NEI", "STA", "STA",
      "NEI", "ENT", "ENT", "ESC", "KLE"))

  # select with two columns
  expect_identical(
    septic_patients[1:10,] %>%
      pull(bactid),
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus, species) %>%
      guess_bactid())
})

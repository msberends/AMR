context("bactid.R")

test_that("as.bactid works", {
  expect_identical(
    as.character(as.bactid(c("E. coli", "H. influenzae"))),
    c("ESCCOL", "HAEINF"))

  expect_equal(as.character(as.bactid("Escherichia coli")), "ESCCOL")
  expect_equal(as.character(as.bactid("P. aer")), "PSEAER") # not Pasteurella aerogenes

  expect_equal(as.character(as.bactid("Negative rods")), "GNR")

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

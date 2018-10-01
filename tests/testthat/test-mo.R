context("mo.R")

test_that("as.mo works", {

  library(dplyr)
  MOs <- AMR::microorganisms %>% filter(!is.na(mo), nchar(mo) > 3)
  expect_identical(as.character(MOs$mo), as.character(as.mo(MOs$mo)))

  expect_identical(
    as.character(as.mo(c("E. coli", "H. influenzae"))),
    c("B_ESCHR_COL", "B_HMPHL_INF"))

  expect_equal(as.character(as.mo("Escherichia coli")), "B_ESCHR_COL")
  expect_equal(as.character(as.mo("Escherichia  coli")), "B_ESCHR_COL")
  expect_equal(as.character(as.mo("Escherichia  species")), "B_ESCHR")
  expect_equal(as.character(as.mo("Escherichia")), "B_ESCHR")
  expect_equal(as.character(as.mo(" B_ESCHR_COL ")), "B_ESCHR_COL")
  expect_equal(as.character(as.mo("e coli")), "B_ESCHR_COL") # not Campylobacter
  expect_equal(as.character(as.mo("klpn")), "B_KLBSL_PNE")
  expect_equal(as.character(as.mo("Klebsiella")), "B_KLBSL")
  expect_equal(as.character(as.mo("K. pneu rhino")), "B_KLBSL_PNE_RHI") # K. pneumoniae subspp. rhinoscleromatis
  expect_equal(as.character(as.mo("Bartonella")), "B_BRTNL")
  expect_equal(as.character(as.mo("C. difficile")), "B_CTRDM_DIF")
  expect_equal(as.character(as.mo("L. pneumophila")), "B_LGNLL_PNE")
  # expect_equal(as.character(as.mo("L. non pneumophila")), "LEGNON")
  # expect_equal(as.character(as.mo("S. beta-haemolytic")), "STCHAE")
  expect_equal(as.character(as.mo("Strepto")), "B_STRPTC")
  expect_equal(as.character(as.mo("Streptococcus")), "B_STRPTC") # not Peptostreptoccus

  expect_equal(as.character(as.mo(c("GAS", "GBS"))), c("B_STRPTC_GRA", "B_STRPTC_GRB"))


  expect_equal(as.character(as.mo("S. pyo")), "B_STRPTC_PYO") # not Actinomyces pyogenes

  expect_equal(as.character(as.mo("P. aer")), "B_PDMNS_AER") # not Pasteurella aerogenes

  # GLIMS
  expect_equal(as.character(as.mo("bctfgr")), "B_BCTRD_FRA")

  expect_equal(as.character(as.mo("MRSE")), "B_STPHY_EPI")
  expect_equal(as.character(as.mo("VRE")), "B_ENTRC")
  expect_equal(as.character(as.mo("MRPA")), "B_PDMNS_AER")
  expect_equal(as.character(as.mo("PISP")), "B_STRPTC_PNE")
  expect_equal(as.character(as.mo("PRSP")), "B_STRPTC_PNE")
  expect_equal(as.character(as.mo("VISP")), "B_STRPTC_PNE")
  expect_equal(as.character(as.mo("VRSP")), "B_STRPTC_PNE")

  expect_equal(as.character(as.mo("CNS")), "B_STPHY_CNS")
  expect_equal(as.character(as.mo("CoNS")), "B_STPHY_CNS")
  expect_equal(as.character(as.mo("CPS")), "B_STPHY_CPS")
  expect_equal(as.character(as.mo("CoPS")), "B_STPHY_CPS")

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
    rep("B_STPHY_AUR", 8))

  # check for Becker classification
  expect_identical(as.character(guess_mo("S. epidermidis", Becker = FALSE)), "B_STPHY_EPI")
  expect_identical(as.character(guess_mo("S. epidermidis", Becker = TRUE)),  "B_STPHY_CNS")
  expect_identical(as.character(guess_mo("STAEPI",         Becker = TRUE)),  "B_STPHY_CNS")
  expect_identical(as.character(guess_mo("S. intermedius", Becker = FALSE)), "B_STRPTC_INT") # Strep (!) intermedius
  expect_identical(as.character(guess_mo("Sta intermedius",Becker = FALSE)), "B_STPHY_INT")
  expect_identical(as.character(guess_mo("Sta intermedius",Becker = TRUE)),  "B_STPHY_CPS")
  expect_identical(as.character(guess_mo("STAINT",         Becker = TRUE)),  "B_STPHY_CPS")
  # aureus must only be influenced if Becker = "all"
  expect_identical(as.character(guess_mo("STAAUR", Becker = FALSE)), "B_STPHY_AUR")
  expect_identical(as.character(guess_mo("STAAUR", Becker = TRUE)),  "B_STPHY_AUR")
  expect_identical(as.character(guess_mo("STAAUR", Becker = "all")), "B_STPHY_CPS")

  # check for Lancefield classification
  expect_identical(as.character(guess_mo("S. pyogenes", Lancefield = FALSE)),    "B_STRPTC_PYO")
  expect_identical(as.character(guess_mo("S. pyogenes", Lancefield = TRUE)),     "B_STRPTC_GRA")
  expect_identical(as.character(guess_mo("STCPYO",      Lancefield = TRUE)),     "B_STRPTC_GRA") # group A
  expect_identical(as.character(guess_mo("S. agalactiae",  Lancefield = FALSE)), "B_STRPTC_AGA")
  expect_identical(as.character(guess_mo("S. agalactiae",  Lancefield = TRUE)),  "B_STRPTC_GRB") # group B
  expect_identical(as.character(guess_mo("S. equisimilis", Lancefield = FALSE)), "B_STRPTC_DYS_EQU")
  expect_identical(as.character(guess_mo("S. equisimilis", Lancefield = TRUE)),  "B_STRPTC_GRC") # group C
  # Enterococci must only be influenced if Lancefield = "all"
  expect_identical(as.character(guess_mo("E. faecium", Lancefield = FALSE)),     "B_ENTRC_IUM")
  expect_identical(as.character(guess_mo("E. faecium", Lancefield = TRUE)),      "B_ENTRC_IUM")
  expect_identical(as.character(guess_mo("E. faecium", Lancefield = "all")),     "B_STRPTC_GRD") # group D
  expect_identical(as.character(guess_mo("S. anginosus",   Lancefield = FALSE)), "B_STRPTC_ANG")
  expect_identical(as.character(guess_mo("S. anginosus",   Lancefield = TRUE)),  "B_STRPTC_GRF") # group F
  expect_identical(as.character(guess_mo("S. sanguinis",   Lancefield = FALSE)), "B_STRPTC_SAN")
  expect_identical(as.character(guess_mo("S. sanguinis",   Lancefield = TRUE)),  "B_STRPTC_GRH") # group H
  expect_identical(as.character(guess_mo("S. salivarius",  Lancefield = FALSE)), "B_STRPTC_SAL")
  expect_identical(as.character(guess_mo("S. salivarius",  Lancefield = TRUE)),  "B_STRPTC_GRK") # group K

  library(dplyr)

  # select with one column
  expect_identical(
    septic_patients[1:10,] %>%
      left_join_microorganisms() %>%
      select(genus) %>%
      as.mo() %>%
      as.character(),
    c("B_ESCHR", "B_ESCHR", "B_STPHY", "B_STPHY", "B_STPHY",
      "B_STPHY", "B_STPHY", "B_STPHY", "B_STPHY", "B_STPHY"))

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
  expect_output(print(as.mo(c("B_ESCHR_COL", NA))))

  # helper function
  expect_identical(as.mo("B_ESCHR_COL"),
                   guess_mo("B_ESCHR_COL"))

  # test pull
  expect_equal(nrow(septic_patients %>% mutate(mo = as.mo(mo))),
               2000)

  # test data.frame
  expect_equal(nrow(data.frame(test = as.mo("B_ESCHR_COL"))),
               1)

  # check empty values
  expect_equal(as.character(suppressWarnings(as.mo(""))),
               NA_character_)

  # check less prevalent MOs
  expect_equal(as.character(as.mo("Gomphosphaeria aponina delicatula")), "B_GMPHS_APO_DEL")
  expect_equal(as.character(as.mo("G apo deli")), "B_GMPHS_APO_DEL")
  expect_equal(as.character(as.mo("Gomphosphaeria  aponina")), "B_GMPHS_APO")
  expect_equal(as.character(as.mo("Gomphosphaeria  species")), "B_GMPHS")
  expect_equal(as.character(as.mo("Gomphosphaeria")), "B_GMPHS")
  expect_equal(as.character(as.mo(" B_GMPHS_APO ")), "B_GMPHS_APO")
  expect_equal(as.character(as.mo("g aponina")), "B_GMPHS_APO")

  # check old names
  expect_equal(suppressMessages(as.character(as.mo("Escherichia blattae"))), "B_SHMWL_BLA")

  # check uncertain names
  expect_equal(suppressWarnings(as.character(as.mo("esco extra_text", allow_uncertain = FALSE))), NA_character_)
  expect_equal(suppressWarnings(as.character(as.mo("esco extra_text", allow_uncertain = TRUE))), "B_ESCHR_COL")
  expect_warning(as.mo("esco extra_text", allow_uncertain = TRUE))

  # predefined reference_df
  expect_equal(as.character(as.mo("TestingOwnID",
                                  reference_df = data.frame(a = "TestingOwnID", b = "B_ESCHR_COL"))),
               "B_ESCHR_COL")
  expect_equal(as.character(as.mo(c("TestingOwnID", "E. coli"),
                                  reference_df = data.frame(a = "TestingOwnID", b = "B_ESCHR_COL"))),
               c("B_ESCHR_COL", "B_ESCHR_COL"))
  expect_warning(as.character(as.mo("TestingOwnID",
                                  reference_df = NULL)))

})

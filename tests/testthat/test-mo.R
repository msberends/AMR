# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

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
  expect_equal(as.character(as.mo(22242416)), "B_ESCHR_COL")
  expect_equal(as.character(as.mo("Escherichia  species")), "B_ESCHR")
  expect_equal(as.character(as.mo("Escherichia")), "B_ESCHR")
  expect_equal(as.character(as.mo("Esch spp.")), "B_ESCHR")
  expect_equal(as.character(as.mo(" B_ESCHR_COL ")), "B_ESCHR_COL")
  expect_equal(as.character(as.mo("e coli")), "B_ESCHR_COL") # not Campylobacter
  expect_equal(as.character(as.mo("klpn")), "B_KLBSL_PNE")
  expect_equal(as.character(as.mo("Klebsiella")), "B_KLBSL")
  expect_equal(as.character(as.mo("K. pneu rhino")), "B_KLBSL_PNE_RHI") # K. pneumoniae subspp. rhinoscleromatis
  expect_equal(as.character(as.mo("Bartonella")), "B_BRTNL")
  expect_equal(as.character(as.mo("C. difficile")), "B_CLSTR_DIF")
  expect_equal(as.character(as.mo("L. pneumophila")), "B_LGNLL_PNE")
  expect_equal(as.character(as.mo("Strepto")), "B_STRPT")
  expect_equal(as.character(as.mo("Streptococcus")), "B_STRPT") # not Peptostreptoccus
  expect_equal(as.character(as.mo("Estreptococos grupo B")), "B_STRPT_GRB")
  expect_equal(as.character(as.mo("Group B Streptococci")), "B_STRPT_GRB")
  expect_equal(as.character(as.mo("B_STRPTC")), "B_STRPT") # old MO code (<=v0.5.0)

  expect_equal(as.character(as.mo(c("GAS", "GBS"))), c("B_STRPT_GRA", "B_STRPT_GRB"))

  expect_equal(as.character(as.mo("S. pyo")), "B_STRPT_PYO") # not Actinomyces pyogenes

  # GLIMS
  expect_equal(as.character(as.mo("bctfgr")), "B_BCTRD_FRA")

  expect_equal(as.character(as.mo("MRSE")), "B_STPHY_EPI")
  expect_equal(as.character(as.mo("VRE")), "B_ENTRC")
  expect_equal(as.character(as.mo("MRPA")), "B_PSDMN_AER")
  expect_equal(as.character(as.mo("PISP")), "B_STRPT_PNE")
  expect_equal(as.character(as.mo("PRSP")), "B_STRPT_PNE")
  expect_equal(as.character(as.mo("VISP")), "B_STRPT_PNE")
  expect_equal(as.character(as.mo("VRSP")), "B_STRPT_PNE")

  expect_equal(as.character(as.mo("CNS")), "B_STPHY_CNS")
  expect_equal(as.character(as.mo("CoNS")), "B_STPHY_CNS")
  expect_equal(as.character(as.mo("CPS")), "B_STPHY_CPS")
  expect_equal(as.character(as.mo("CoPS")), "B_STPHY_CPS")
  expect_equal(as.character(as.mo("VGS")), "B_STRPT_VIR")
  expect_equal(as.character(as.mo("streptococcus milleri")), "B_STRPT_MIL")
  

  expect_equal(as.character(as.mo(c("Gram negative", "Gram positive"))), c("B_GRAMN", "B_GRAMP"))

  # prevalent MO
  expect_identical(
    suppressWarnings(as.character(
      as.mo(c("stau",
              "STAU",
              "staaur",
              "S. aureus",
              "S aureus",
              "Sthafilokkockus aureeuzz",
              "Staphylococcus aureus",
              "MRSA",
              "VISA")))),
    rep("B_STPHY_AUR", 9))
  expect_identical(
    as.character(
      as.mo(c('EHEC', 'EPEC', 'EIEC', 'STEC', 'ATEC', 'UPEC'))),
    rep("B_ESCHR_COL", 6))
  # unprevalent MO
  expect_identical(
    as.character(
      as.mo(c("burnod",
              "B. nodosa",
              "B nodosa",
              "Burkholderia nodosa"))),
    rep("B_BRKHL_NOD", 4))

  # empty values
  expect_identical(as.character(as.mo(c("", NA, NaN))), rep(NA_character_, 3))
  # too few characters
  expect_warning(as.mo("ab"))

  expect_equal(suppressWarnings(as.character(as.mo(c("Qq species", "", "CRS", "K. pneu rhino", "esco")))),
               c("UNKNOWN", NA_character_, "B_STNTR_MAL", "B_KLBSL_PNE_RHI", "B_ESCHR_COL"))

  # check for Becker classification
  expect_identical(as.character(as.mo("S. epidermidis", Becker = FALSE)), "B_STPHY_EPI")
  expect_identical(as.character(as.mo("S. epidermidis", Becker = TRUE)),  "B_STPHY_CNS")
  expect_identical(as.character(as.mo("STAEPI",         Becker = TRUE)),  "B_STPHY_CNS")
  expect_identical(as.character(as.mo("S. intermedius", Becker = FALSE)), "B_STPHY_INT")
  expect_identical(as.character(as.mo("Sta intermedius",Becker = FALSE)), "B_STPHY_INT")
  expect_identical(as.character(as.mo("Sta intermedius",Becker = TRUE)),  "B_STPHY_CPS")
  expect_identical(as.character(as.mo("STAINT",         Becker = TRUE)),  "B_STPHY_CPS")
  # aureus must only be influenced if Becker = "all"
  expect_identical(as.character(as.mo("STAAUR", Becker = FALSE)), "B_STPHY_AUR")
  expect_identical(as.character(as.mo("STAAUR", Becker = TRUE)),  "B_STPHY_AUR")
  expect_identical(as.character(as.mo("STAAUR", Becker = "all")), "B_STPHY_CPS")

  # check for Lancefield classification
  expect_identical(as.character(as.mo("S. pyogenes", Lancefield = FALSE)),    "B_STRPT_PYO")
  expect_identical(as.character(as.mo("S. pyogenes", Lancefield = TRUE)),     "B_STRPT_GRA")
  expect_identical(as.character(as.mo("STCPYO",      Lancefield = TRUE)),     "B_STRPT_GRA") # group A
  expect_identical(as.character(as.mo("S. agalactiae",  Lancefield = FALSE)), "B_STRPT_AGA")
  expect_identical(as.character(as.mo("S. agalactiae",  Lancefield = TRUE)),  "B_STRPT_GRB") # group B
  expect_identical(as.character(suppressWarnings(as.mo("estreptococos grupo B"))), "B_STRPT_GRB")
  expect_identical(as.character(as.mo("S. equisimilis", Lancefield = FALSE)), "B_STRPT_DYS_EQU")
  expect_identical(as.character(as.mo("S. equisimilis", Lancefield = TRUE)),  "B_STRPT_GRC") # group C
  # Enterococci must only be influenced if Lancefield = "all"
  expect_identical(as.character(as.mo("E. faecium", Lancefield = FALSE)),     "B_ENTRC_IUM")
  expect_identical(as.character(as.mo("E. faecium", Lancefield = TRUE)),      "B_ENTRC_IUM")
  expect_identical(as.character(as.mo("E. faecium", Lancefield = "all")),     "B_STRPT_GRD") # group D
  expect_identical(as.character(as.mo("S. anginosus",   Lancefield = FALSE)), "B_STRPT_ANG")
  expect_identical(as.character(as.mo("S. anginosus",   Lancefield = TRUE)),  "B_STRPT_GRF") # group F
  expect_identical(as.character(as.mo("S. sanguinis",   Lancefield = FALSE)), "B_STRPT_SAN")
  expect_identical(as.character(as.mo("S. sanguinis",   Lancefield = TRUE)),  "B_STRPT_GRH") # group H
  expect_identical(as.character(as.mo("S. salivarius",  Lancefield = FALSE)), "B_STRPT_SAL")
  expect_identical(as.character(as.mo("S. salivarius",  Lancefield = TRUE)),  "B_STRPT_GRK") # group K

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
      as.mo())

  # unknown results
  expect_warning(as.mo(c("INVALID", "Yeah, unknown")))

  # too many columns
  expect_error(septic_patients %>% select(1:3) %>% as.mo())

  # print
  expect_output(print(as.mo(c("B_ESCHR_COL", NA))))

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
  expect_equal(as.character(as.mo("Gomphosphaeria apo del")), "B_GMPHS_APO_DEL")
  expect_equal(as.character(as.mo("G apo deli")), "B_GMPHS_APO_DEL")
  expect_equal(as.character(as.mo("Gomphosphaeria  aponina")), "B_GMPHS_APO")
  expect_equal(as.character(as.mo("Gomphosphaeria  species")), "B_GMPHS")
  expect_equal(as.character(as.mo("Gomphosphaeria")), "B_GMPHS")
  expect_equal(as.character(as.mo(" B_GMPHS_APO ")), "B_GMPHS_APO")
  expect_equal(as.character(as.mo("g aponina")), "B_GMPHS_APO")

  # check old names
  expect_equal(suppressMessages(as.character(as.mo("Escherichia blattae"))), "B_SHMWL_BLA")
  print(mo_renamed())

  # check uncertain names
  expect_equal(suppressWarnings(as.character(as.mo("staaur extratest", allow_uncertain = TRUE))), "B_STPHY_AUR")
  expect_equal(suppressWarnings(as.character(as.mo("staaur extratest", allow_uncertain = FALSE))), "UNKNOWN")
  expect_warning(as.mo("esco extra_text", allow_uncertain = TRUE))
  expect_equal(suppressWarnings(as.character(as.mo("unexisting aureus", allow_uncertain = 3))), "B_STPHY_AUR")
  expect_equal(suppressWarnings(as.character(as.mo("unexisting staphy", allow_uncertain = 3))), "B_STPHY")
  expect_equal(suppressWarnings(as.character(as.mo("Staphylococcus aureus unexisting", allow_uncertain = 3))), "B_STPHY_AUR")

  # predefined reference_df
  expect_equal(as.character(as.mo("TestingOwnID",
                                  reference_df = data.frame(mycol = "TestingOwnID", mo = "B_ESCHR_COL"))),
               "B_ESCHR_COL")
  expect_equal(as.character(as.mo(c("TestingOwnID", "E. coli"),
                                  reference_df = data.frame(mycol = "TestingOwnID", mo = "B_ESCHR_COL"))),
               c("B_ESCHR_COL", "B_ESCHR_COL"))
  expect_warning(as.mo("TestingOwnID", reference_df = NULL))
  expect_error(as.mo("E. coli", reference_df = data.frame(mycol = "TestingOwnID")))

  # combination of existing mo and other code
  expect_identical(as.character(as.mo(c("B_ESCHR_COL", "ESCCOL"))),
                   c("B_ESCHR_COL", "B_ESCHR_COL"))

  # expect_equal(mo_fullname(c("E. spp.",
  #                            "E. spp",
  #                            "E. species")),
  #              rep("Escherichia species", 3))

  # from different sources
  expect_equal(as.character(as.mo(
    c("PRTMIR", "bclcer", "B_ESCHR_COL"))),
    c("B_PROTS_MIR", "B_BCLLS_CER", "B_ESCHR_COL"))

  # hard to find
  expect_equal(as.character(suppressWarnings(as.mo(
    c("Microbacterium paraoxidans",
      "Streptococcus suis (bovis gr)",
      "Raoultella (here some text) terrigena")))),
    c("B_MCRBC_PAR", "B_STRPT_SUI", "B_RLTLL_TER"))
  print(mo_uncertainties())

  # Salmonella (City) are all actually Salmonella enterica spp (City)
  expect_equal(as.character(suppressWarnings(as.mo("Salmonella Goettingen"))),
               "B_SLMNL_ENT")
  expect_equal(as.character(as.mo("Salmonella Group A")), "B_SLMNL")

  # no virusses
  expect_warning(as.mo("Virus"))

  # summary
  expect_equal(length(summary(septic_patients$mo)), 6)

  # WHONET codes and NA/NaN
  expect_equal(as.character(as.mo(c("xxx", "na", "nan"), debug = TRUE)),
               rep(NA_character_, 3))
  expect_equal(as.character(as.mo("con")), "UNKNOWN")
  expect_equal(as.character(as.mo("xxx")), NA_character_)
  expect_equal(as.character(as.mo(c("xxx", "con", "eco"))), c(NA_character_, "UNKNOWN", "B_ESCHR_COL"))
    expect_equal(as.character(as.mo(c("other", "none", "unknown"))),
               rep("UNKNOWN", 3))

  expect_null(mo_failures())
  expect_true(septic_patients %>% pull(mo) %>% is.mo())

  # expect_equal(get_mo_code("test", "mo"), "test")
  # expect_equal(length(get_mo_code("Escherichia", "genus")),
  #              nrow(AMR::microorganisms[base::which(AMR::microorganisms[, "genus"] %in% "Escherichia"),]))

  expect_error(translate_allow_uncertain(5))

  # very old MO codes (<= v0.5.0)
  expect_equal(as.character(as.mo("F_CCCCS_NEO")), "F_CRYPT_NEO")
  expect_equal(as.character(as.mo("F_CANDD_GLB")), "F_CANDD_GLA")
  
  # debug mode
  expect_output(print(suppressMessages(suppressWarnings(as.mo("kshgcjkhsdgkshjdfsfvsdfv", debug = TRUE, allow_uncertain = 3)))))

  # ..coccus
  expect_equal(as.character(as.mo(c("meningococ", "gonococ", "pneumococ"))), 
               c("B_NESSR_MEN", "B_NESSR_GON", "B_STRPT_PNE"))
  # yeasts and fungi
  expect_equal(suppressWarnings(as.character(as.mo(c("yeasts", "fungi")))), 
               c("F_YEAST", "F_FUNGUS"))
  
  # print tibble
  expect_output(print(tibble(mo = as.mo("B_STRPT_PNE"))))
  
  # assigning and subsetting
  x <- septic_patients$mo
  expect_s3_class(x[1], "mo")
  expect_s3_class(x[[1]], "mo")
  expect_s3_class(c(x[1], x[9]), "mo")
  expect_warning(x[1] <- "invalid code")
  expect_warning(x[[1]] <- "invalid code")
  expect_warning(c(x[1], "test"))
})

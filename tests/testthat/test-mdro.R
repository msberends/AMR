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

context("mdro.R")

test_that("mdro works", {
  
  skip_on_cran()
  
  library(dplyr)

  expect_error(suppressWarnings(mdro(example_isolates, country = "invalid", col_mo = "mo", info = TRUE)))
  expect_error(suppressWarnings(mdro(example_isolates, country = "fr", info = TRUE)))
  expect_error(mdro(example_isolates, guideline = c("BRMO", "MRGN"), info = TRUE))
  expect_error(mdro(example_isolates, col_mo = "invalid", info = TRUE))

  outcome <- mdro(example_isolates)
  outcome <- eucast_exceptional_phenotypes(example_isolates, info = TRUE)
  # check class
  expect_equal(outcome %>% class(), c("ordered", "factor"))

  outcome <- mdro(example_isolates, "nl", info = TRUE)
  # check class
  expect_equal(outcome %>% class(), c("ordered", "factor"))

  # example_isolates should have these finding using Dutch guidelines
  expect_equal(outcome %>% freq() %>% pull(count),
               c(1969, 25, 6)) # 1969 neg, 25 unconfirmed, 6 pos

  expect_equal(brmo(example_isolates, info = FALSE),
               mdro(example_isolates, guideline = "BRMO", info = FALSE))

  # test Dutch P. aeruginosa MDRO
  expect_equal(
    as.character(mdro(data.frame(mo = as.mo("P. aeruginosa"),
                                 cfta = "S",
                                 cipr = "S",
                                 mero = "S",
                                 imip = "S",
                                 gent = "S",
                                 tobr = "S",
                                 pita = "S"),
                      guideline = "BRMO",
                      col_mo = "mo",
                      info = FALSE)),
    "Negative")
  expect_equal(
    as.character(mdro(data.frame(mo = as.mo("P. aeruginosa"),
                                 cefta = "R",
                                 cipr = "R",
                                 mero = "R",
                                 imip = "R",
                                 gent = "R",
                                 tobr = "R",
                                 pita = "R"),
                      guideline = "BRMO",
                      col_mo = "mo",
                      info = FALSE)),
    "Positive")
  
  # German 3MRGN and 4MRGN
  expect_equal(as.character(mrgn(
    data.frame(mo = c("E. coli", "E. coli", "K. pneumoniae", "E. coli",
                      "A. baumannii", "A. baumannii", "A. baumannii",
                      "P. aeruginosa", "P. aeruginosa", "P. aeruginosa"), 
               PIP = c("S", "R", "R", "S",
                       "S", "R", "R",
                       "S", "R", "R"),
               CTX = c("S", "R", "R", "S",
                       "R", "R", "R",
                       "R", "R", "R"),
               IPM = c("S", "R", "S", "R",
                       "R", "R", "S",
                       "S", "R", "R"),
               CIP = c("S", "R", "R", "S",
                       "R", "R", "R",
                       "R", "S", "R"),
               stringsAsFactors = FALSE))),
    c("Negative", "4MRGN", "3MRGN", "4MRGN",  "4MRGN",  "4MRGN", "3MRGN", "Negative", "3MRGN", "4MRGN"))
  
  # MDR TB
  expect_equal(
    # select only rifampicine, mo will be determined automatically (as M. tuberculosis),
    # number of mono-resistant strains should be equal to number of rifampicine-resistant strains
    example_isolates %>% select(RIF) %>% mdr_tb() %>% freq() %>% pull(count) %>% .[2],
    count_R(example_isolates$RIF))

  sample_rsi <- function() {
    sample(c("S", "I", "R"),
           size = 5000,
           prob = c(0.5, 0.1, 0.4),
           replace = TRUE)
  }
  expect_gt(
    #suppressWarnings(
      data.frame(rifampicin = sample_rsi(),
                 inh = sample_rsi(),
                 gatifloxacin = sample_rsi(),
                 eth = sample_rsi(),
                 pza = sample_rsi(),
                 MFX = sample_rsi(),
                 KAN = sample_rsi()) %>%
        mdr_tb() %>%
        n_distinct()
      #)
      ,
    2)
  
  # check the guideline by Magiorakos  et al. (2012), the default guideline
  stau <- tribble(
    ~mo,         ~GEN, ~RIF, ~CPT, ~OXA, ~CIP, ~MFX, ~SXT, ~FUS, ~VAN, ~TEC, ~TLV, ~TGC, ~CLI, ~DAP, ~ERY, ~LNZ, ~CHL, ~FOS, ~QDA, ~TCY, ~DOX, ~MNO,
    "S. aureus",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "S. aureus",  "R",  "R",  "R",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "S. aureus",  "S",  "S",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",
    "S. aureus",  "R",  "R",  "I",  "I",  "I",  "I",  "I",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R"
  )
  expect_equal(as.integer(mdro(stau)), c(1:4))
  expect_s3_class(mdro(stau, verbose = TRUE), "data.frame")
  
  ente <- tribble(
    ~mo,            ~GEH, ~STH, ~IPM, ~MEM, ~DOR, ~CIP, ~LVX, ~MFX, ~VAN, ~TEC, ~TGC, ~DAP, ~LNZ, ~AMP, ~QDA, ~DOX, ~MNO,
    "Enterococcus",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "Enterococcus",  "R",  "R",  "R",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "Enterococcus",  "S",  "S",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",
    "Enterococcus",  "R",  "R",  "I",  "I",  "I",  "I",  "I",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R"
  )
  expect_equal(as.integer(mdro(ente)), c(1:4))
  expect_s3_class(mdro(ente, verbose = TRUE), "data.frame")
  
  entero <- tribble(
    ~mo,       ~GEN, ~TOB, ~AMK, ~NET, ~CPT, ~TCC, ~TZP, ~ETP, ~IPM, ~MEM, ~DOR, ~CZO, ~CXM, ~CTX, ~CAZ, ~FEP, ~FOX, ~CTT, ~CIP, ~SXT, ~TGC, ~ATM, ~AMP, ~AMC, ~SAM, ~CHL, ~FOS, ~COL, ~TCY, ~DOX, ~MNO,
    "E. coli",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "E. coli",  "R",  "R",  "R",  "R",  "R",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "E. coli",  "S",  "S",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",
    "E. coli",  "R",  "R",  "I",  "I",  "I",  "I",  "I",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R"
  )
  expect_equal(as.integer(mdro(entero)), c(1:4))
  expect_s3_class(mdro(entero, verbose = TRUE), "data.frame")
  
  pseud <- tribble(
    ~mo,             ~GEN, ~TOB, ~AMK, ~NET, ~IPM, ~MEM, ~DOR, ~CAZ, ~FEP, ~CIP, ~LVX, ~TCC, ~TZP, ~ATM, ~FOS, ~COL, ~PLB,
    "P. aeruginosa",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "P. aeruginosa",  "R",  "S",  "S",  "S",  "R",  "S",  "S",  "S",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "P. aeruginosa",  "S",  "S",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",
    "P. aeruginosa",  "R",  "R",  "I",  "I",  "I",  "I",  "I",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R"
  )
  expect_equal(as.integer(mdro(pseud)), c(1:4))
  expect_s3_class(mdro(pseud, verbose = TRUE), "data.frame")
  
  acin <- tribble(
    ~mo,            ~GEN, ~TOB, ~AMK, ~NET, ~IPM, ~MEM, ~DOR, ~CIP, ~LVX, ~TZP, ~TCC, ~CTX, ~CRO, ~CAZ, ~FEP, ~SXT, ~SAM, ~COL, ~PLB, ~TCY, ~DOX, ~MNO,
    "A. baumannii",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "A. baumannii",  "R",  "R",  "R",  "R",  "S",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "S",  "R",  "S",  "S",  "S",  "S",  "S",  "S",  "S",
    "A. baumannii",  "S",  "S",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",
    "A. baumannii",  "R",  "R",  "I",  "I",  "I",  "I",  "I",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R",  "R"
  )
  expect_equal(as.integer(mdro(acin)), c(1:4))
  expect_s3_class(mdro(acin, verbose = TRUE), "data.frame")
  
})

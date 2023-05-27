# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

expect_error(mdro(example_isolates, guideline = c("BRMO", "MRGN"), info = TRUE))
expect_error(mdro(example_isolates, col_mo = "invalid", info = TRUE))

expect_stdout(suppressMessages(suppressWarnings(mdro(example_isolates, info = TRUE))))
expect_stdout(suppressMessages(suppressWarnings(mdro(example_isolates, "eucast3.1", info = TRUE))))
expect_stdout(suppressMessages(suppressWarnings(mdro(example_isolates, "eucast3.2", info = TRUE))))
expect_stdout(suppressMessages(suppressWarnings(mdro(example_isolates, "eucast3.3", info = TRUE))))
expect_stdout(outcome <- suppressMessages(suppressWarnings(eucast_exceptional_phenotypes(example_isolates, info = TRUE))))
# check class
expect_identical(class(outcome), c("ordered", "factor"))

expect_stdout(outcome <- mdro(example_isolates, "nl", info = TRUE))
# check class
expect_identical(class(outcome), c("ordered", "factor"))

# example_isolates should have these finding using Dutch guidelines
expect_equal(
  as.double(table(outcome)),
  c(1954, 24, 6)
) # 1954 neg, 24 unconfirmed, 6 pos, rest is NA

expect_equal(
  brmo(example_isolates, info = FALSE),
  mdro(example_isolates, guideline = "BRMO", info = FALSE)
)

# test Dutch P. aeruginosa MDRO
expect_equal(
  as.character(mdro(data.frame(
    mo = as.mo("P. aeruginosa"),
    cfta = "S",
    cipr = "S",
    mero = "S",
    imip = "S",
    gent = "S",
    tobr = "S",
    pita = "S"
  ),
  guideline = "BRMO",
  col_mo = "mo",
  info = FALSE
  )),
  "Negative"
)
expect_equal(
  as.character(mdro(data.frame(
    mo = as.mo("P. aeruginosa"),
    cefta = "R",
    cipr = "R",
    mero = "R",
    imip = "R",
    gent = "R",
    tobr = "R",
    pita = "R"
  ),
  guideline = "BRMO",
  col_mo = "mo",
  info = FALSE
  )),
  "Positive"
)

# German 3MRGN and 4MRGN
expect_equal(
  as.character(mrgn(
    data.frame(
      mo = c(
        "E. coli", "E. coli", "K. pneumoniae", "E. coli",
        "A. baumannii", "A. baumannii", "A. baumannii",
        "P. aeruginosa", "P. aeruginosa", "P. aeruginosa"
      ),
      PIP = c(
        "S", "R", "R", "S",
        "S", "R", "R",
        "S", "R", "R"
      ),
      CTX = c(
        "S", "R", "R", "S",
        "R", "R", "R",
        "R", "R", "R"
      ),
      IPM = c(
        "S", "R", "S", "R",
        "R", "R", "S",
        "S", "R", "R"
      ),
      CIP = c(
        "S", "R", "R", "S",
        "R", "R", "R",
        "R", "S", "R"
      ),
      stringsAsFactors = FALSE
    )
  )),
  c("Negative", "4MRGN", "3MRGN", "4MRGN", "4MRGN", "4MRGN", "3MRGN", "Negative", "3MRGN", "4MRGN")
)

# MDR TB
expect_equal(
  # select only rifampicine, mo will be determined automatically (as M. tuberculosis),
  # number of mono-resistant strains should be equal to number of rifampicine-resistant strains
  as.double(table(mdr_tb(example_isolates[, "RIF", drop = FALSE])))[2],
  count_R(example_isolates$RIF)
)

x <- data.frame(
  rifampicin = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  inh = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  gatifloxacin = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  eth = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  pza = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  MFX = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5)),
  KAN = random_sir(5000, prob_sir = c(0.4, 0.1, 0.5))
)
expect_true(length(unique(mdr_tb(x))) > 2)

# check the guideline by Magiorakos  et al. (2012), the default guideline
stau <- data.frame(
  mo = c("S. aureus", "S. aureus", "S. aureus", "S. aureus"),
  GEN = c("R", "R", "S", "R"),
  RIF = c("S", "R", "S", "R"),
  CPT = c("S", "R", "R", "R"),
  OXA = c("S", "R", "R", "R"),
  CIP = c("S", "S", "R", "R"),
  MFX = c("S", "S", "R", "R"),
  SXT = c("S", "S", "R", "R"),
  FUS = c("S", "S", "R", "R"),
  VAN = c("S", "S", "R", "R"),
  TEC = c("S", "S", "R", "R"),
  TLV = c("S", "S", "R", "R"),
  TGC = c("S", "S", "R", "R"),
  CLI = c("S", "S", "R", "R"),
  DAP = c("S", "S", "R", "R"),
  ERY = c("S", "S", "R", "R"),
  LNZ = c("S", "S", "R", "R"),
  CHL = c("S", "S", "R", "R"),
  FOS = c("S", "S", "R", "R"),
  QDA = c("S", "S", "R", "R"),
  TCY = c("S", "S", "R", "R"),
  DOX = c("S", "S", "R", "R"),
  MNO = c("S", "S", "R", "R"),
  stringsAsFactors = FALSE
)
expect_equal(as.integer(mdro(stau)), c(1:4))
expect_inherits(mdro(stau, verbose = TRUE), "data.frame")

ente <- data.frame(
  mo = c("Enterococcus", "Enterococcus", "Enterococcus", "Enterococcus"),
  GEH = c("R", "R", "S", "R"),
  STH = c("S", "R", "S", "R"),
  IPM = c("S", "R", "R", "R"),
  MEM = c("S", "R", "R", "R"),
  DOR = c("S", "S", "R", "R"),
  CIP = c("S", "S", "R", "R"),
  LVX = c("S", "S", "R", "R"),
  MFX = c("S", "S", "R", "R"),
  VAN = c("S", "S", "R", "R"),
  TEC = c("S", "S", "R", "R"),
  TGC = c("S", "S", "R", "R"),
  DAP = c("S", "S", "R", "R"),
  LNZ = c("S", "S", "R", "R"),
  AMP = c("S", "S", "R", "R"),
  QDA = c("S", "S", "R", "R"),
  DOX = c("S", "S", "R", "R"),
  MNO = c("S", "S", "R", "R"),
  stringsAsFactors = FALSE
)
expect_equal(as.integer(mdro(ente)), c(1:4))
expect_inherits(mdro(ente, verbose = TRUE), "data.frame")

entero <- data.frame(
  mo = c("E. coli", "E. coli", "E. coli", "E. coli"),
  GEN = c("R", "R", "S", "R"), TOB = c("S", "R", "S", "R"),
  AMK = c("S", "R", "R", "R"), NET = c("S", "R", "R", "R"),
  CPT = c("S", "R", "R", "R"), TCC = c("S", "R", "R", "R"),
  TZP = c("S", "S", "R", "R"), ETP = c("S", "S", "R", "R"),
  IPM = c("S", "S", "R", "R"), MEM = c("S", "S", "R", "R"),
  DOR = c("S", "S", "R", "R"), CZO = c("S", "S", "R", "R"),
  CXM = c("S", "S", "R", "R"), CTX = c("S", "S", "R", "R"),
  CAZ = c("S", "S", "R", "R"), FEP = c("S", "S", "R", "R"),
  FOX = c("S", "S", "R", "R"), CTT = c("S", "S", "R", "R"),
  CIP = c("S", "S", "R", "R"), SXT = c("S", "S", "R", "R"),
  TGC = c("S", "S", "R", "R"), ATM = c("S", "S", "R", "R"),
  AMP = c("S", "S", "R", "R"), AMC = c("S", "S", "R", "R"),
  SAM = c("S", "S", "R", "R"), CHL = c("S", "S", "R", "R"),
  FOS = c("S", "S", "R", "R"), COL = c("S", "S", "R", "R"),
  TCY = c("S", "S", "R", "R"), DOX = c("S", "S", "R", "R"),
  MNO = c("S", "S", "R", "R"),
  stringsAsFactors = FALSE
)
expect_equal(as.integer(mdro(entero)), c(1:4))
expect_inherits(mdro(entero, verbose = TRUE), "data.frame")

pseud <- data.frame(
  mo = c("P. aeruginosa", "P. aeruginosa", "P. aeruginosa", "P. aeruginosa"),
  GEN = c("R", "R", "S", "R"), TOB = c("S", "S", "S", "R"),
  AMK = c("S", "S", "R", "R"), NET = c("S", "S", "R", "R"),
  IPM = c("S", "R", "R", "R"), MEM = c("S", "S", "R", "R"),
  DOR = c("S", "S", "R", "R"), CAZ = c("S", "S", "R", "R"),
  FEP = c("S", "R", "R", "R"), CIP = c("S", "S", "R", "R"),
  LVX = c("S", "S", "R", "R"), TCC = c("S", "S", "R", "R"),
  TZP = c("S", "S", "R", "R"), ATM = c("S", "S", "R", "R"),
  FOS = c("S", "S", "R", "R"), COL = c("S", "S", "R", "R"),
  PLB = c("S", "S", "R", "R"),
  stringsAsFactors = FALSE
)
expect_equal(as.integer(mdro(pseud)), c(1:4))
expect_inherits(mdro(pseud, verbose = TRUE), "data.frame")

acin <- data.frame(
  mo = c("A. baumannii", "A. baumannii", "A. baumannii", "A. baumannii"),
  GEN = c("R", "R", "S", "R"), TOB = c("S", "R", "S", "R"),
  AMK = c("S", "R", "R", "R"), NET = c("S", "R", "R", "R"),
  IPM = c("S", "S", "R", "R"), MEM = c("S", "R", "R", "R"),
  DOR = c("S", "S", "R", "R"), CIP = c("S", "S", "R", "R"),
  LVX = c("S", "S", "R", "R"), TZP = c("S", "S", "R", "R"),
  TCC = c("S", "S", "R", "R"), CTX = c("S", "S", "R", "R"),
  CRO = c("S", "S", "R", "R"), CAZ = c("S", "S", "R", "R"),
  FEP = c("S", "R", "R", "R"), SXT = c("S", "S", "R", "R"),
  SAM = c("S", "S", "R", "R"), COL = c("S", "S", "R", "R"),
  PLB = c("S", "S", "R", "R"), TCY = c("S", "S", "R", "R"),
  DOX = c("S", "S", "R", "R"), MNO = c("S", "S", "R", "R"),
  stringsAsFactors = FALSE
)
expect_equal(as.integer(mdro(acin)), c(1:4))
expect_inherits(mdro(acin, verbose = TRUE), "data.frame")

# custom rules
custom <- custom_mdro_guideline("CIP == 'R' & age > 60" ~ "Elderly Type A",
  "ERY == 'R' & age > 60" ~ "Elderly Type B",
  as_factor = TRUE
)
expect_stdout(print(custom))
expect_stdout(print(c(custom, custom)))
expect_stdout(print(as.list(custom, custom)))

expect_stdout(x <- mdro(example_isolates, guideline = custom, info = TRUE))
expect_equal(as.double(table(x)), c(1070, 198, 732))

expect_stdout(print(custom_mdro_guideline(AMX == "R" ~ "test", as_factor = FALSE)))
expect_error(custom_mdro_guideline())
expect_error(custom_mdro_guideline("test"))
expect_error(custom_mdro_guideline("test" ~ c(1:3)))
expect_error(custom_mdro_guideline("test" ~ A))
# expect_warning(mdro(example_isolates,
#   # since `test` gives an error, it will be ignored with a warning
#   guideline = custom_mdro_guideline(test ~ "A"),
#   info = FALSE
# ))

# print groups
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
  expect_stdout(x <- mdro(example_isolates %>% group_by(ward), info = TRUE, pct_required_classes = 0))
  expect_stdout(x <- mdro(example_isolates %>% group_by(ward), guideline = custom, info = TRUE))
}

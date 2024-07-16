# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
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

# thoroughly check input table
expect_equal(
  colnames(AMR:::EUCAST_RULES_DF),
  c(
    "if_mo_property", "like.is.one_of", "this_value",
    "and_these_antibiotics", "have_these_values",
    "then_change_these_antibiotics", "to_value",
    "reference.rule", "reference.rule_group",
    "reference.version",
    "note"
  )
)
MOs_mentioned <- unique(AMR:::EUCAST_RULES_DF$this_value)
MOs_mentioned <- sort(trimws(unlist(strsplit(MOs_mentioned[!AMR:::is_valid_regex(MOs_mentioned)], ",", fixed = TRUE))))
MOs_test <- suppressWarnings(suppressMessages(mo_name(MOs_mentioned, keep_synonyms = TRUE, language = NULL)))
expect_true(length(MOs_mentioned[MOs_test != MOs_mentioned]) == 0)

expect_error(suppressWarnings(eucast_rules(example_isolates, col_mo = "Non-existing")))
expect_error(eucast_rules(x = "text"))
expect_error(eucast_rules(data.frame(a = "test")))
expect_error(eucast_rules(data.frame(mo = "test"), rules = "invalid rules set"))

# expect_warning(eucast_rules(data.frame(mo = "Escherichia coli", vancomycin = "S", stringsAsFactors = TRUE)))

expect_identical(
  colnames(example_isolates),
  colnames(suppressWarnings(eucast_rules(example_isolates, info = FALSE)))
)

expect_stdout(suppressMessages(eucast_rules(example_isolates, info = TRUE)))

a <- data.frame(
  mo = c(
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa",
    "Enterobacter cloacae"
  ),
  amox = "-", # Amoxicillin
  stringsAsFactors = FALSE
)
b <- data.frame(
  mo = c(
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa",
    "Enterobacter cloacae"
  ),
  amox = "R", # Amoxicillin
  stringsAsFactors = FALSE
)
expect_identical(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)
expect_stdout(suppressMessages(suppressWarnings(eucast_rules(a, "mo", info = TRUE))))

a <- data.frame(
  mo = c(
    "Staphylococcus aureus",
    "Streptococcus group A"
  ),
  COL = "-", # Colistin
  stringsAsFactors = FALSE
)
b <- data.frame(
  mo = c(
    "Staphylococcus aureus",
    "Streptococcus group A"
  ),
  COL = "R", # Colistin
  stringsAsFactors = FALSE
)
expect_equal(suppressWarnings(eucast_rules(a, "mo", info = FALSE)), b)

# piperacillin must be R in Enterobacteriaceae when tica is R
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
  expect_equal(
    suppressWarnings(
      example_isolates %>%
        filter(mo_family(mo) == "Enterobacteriaceae") %>%
        mutate(
          TIC = as.sir("R"),
          PIP = as.sir("S")
        ) %>%
        eucast_rules(col_mo = "mo", version_expertrules = 3.1, info = FALSE) %>%
        pull(PIP) %>%
        unique() %>%
        as.character()
    ),
    "R"
  )
}

# azithromycin and clarythromycin must be equal to Erythromycin
a <- suppressWarnings(as.sir(eucast_rules(data.frame(
  mo = example_isolates$mo,
  ERY = example_isolates$ERY,
  AZM = as.sir("R"),
  CLR = factor("R"),
  stringsAsFactors = FALSE
),
version_expertrules = 3.1,
only_sir_columns = FALSE
)$CLR))
b <- example_isolates$ERY
expect_identical(
  a[!is.na(b)],
  b[!is.na(b)]
)

# amox is inferred by benzylpenicillin in Kingella kingae
expect_equal(
  suppressWarnings(
    as.list(eucast_rules(
      data.frame(
        mo = as.mo("Kingella kingae"),
        PEN = "S",
        AMX = "-",
        stringsAsFactors = FALSE
      ),
      info = FALSE
    ))$AMX
  ),
  "S"
)

# also test norf
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
  expect_stdout(suppressWarnings(eucast_rules(example_isolates %>% mutate(NOR = "S", NAL = "S"), info = TRUE)))
}

# check verbose output
expect_stdout(suppressWarnings(eucast_rules(example_isolates, verbose = TRUE, rules = "all", info = TRUE)))

# AmpC de-repressed cephalo mutants

expect_identical(
  eucast_rules(data.frame(
    mo = c("Escherichia coli", "Enterobacter cloacae"),
    cefotax = as.sir(c("S", "S"))
  ),
  ampc_cephalosporin_resistance = TRUE,
  info = FALSE
  )$cefotax,
  as.sir(c("S", "R"))
)

expect_identical(
  eucast_rules(data.frame(
    mo = c("Escherichia coli", "Enterobacter cloacae"),
    cefotax = as.sir(c("S", "S"))
  ),
  ampc_cephalosporin_resistance = NA,
  info = FALSE
  )$cefotax,
  as.sir(c("S", NA))
)

expect_identical(
  eucast_rules(data.frame(
    mo = c("Escherichia coli", "Enterobacter cloacae"),
    cefotax = as.sir(c("S", "S"))
  ),
  ampc_cephalosporin_resistance = NULL,
  info = FALSE
  )$cefotax,
  as.sir(c("S", "S"))
)

# EUCAST dosage -----------------------------------------------------------
expect_equal(nrow(eucast_dosage(c("tobra", "genta", "cipro"))), 3)
expect_inherits(eucast_dosage(c("tobra", "genta", "cipro")), "data.frame")



x <- custom_eucast_rules(
  AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
  AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I",
  AMX == "S" ~ AMC == "S"
)
expect_stdout(print(x))
expect_stdout(print(c(x, x)))
expect_stdout(print(as.list(x, x)))

# this custom rules makes 8 changes
expect_equal(nrow(eucast_rules(example_isolates,
  rules = "custom",
  custom_rules = x,
  info = FALSE,
  verbose = TRUE
)),
8,
tolerance = 0.5
)

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
# doi:10.18637/jss.v104.i03                                            #
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

# we must only have EUCAST and CLSI, because otherwise the rules in as.sir() will fail
expect_identical(
  unique(gsub("[^A-Z]", "", AMR::clinical_breakpoints$guideline)),
  c("EUCAST", "CLSI")
)

expect_true(as.sir("S") < as.sir("I"))
expect_true(as.sir("I") < as.sir("R"))
expect_true(is.sir(as.sir("S")))
x <- example_isolates$AMX
expect_inherits(x[1], "sir")
expect_inherits(x[[1]], "sir")
expect_inherits(c(x[1], x[9]), "sir")
expect_inherits(unique(x[1], x[9]), "sir")
pdf(NULL) # prevent Rplots.pdf being created
expect_silent(barplot(as.sir(c("S", "I", "R"))))
expect_silent(plot(as.sir(c("S", "I", "R"))))
if (AMR:::pkg_is_available("ggplot2")) {
  expect_inherits(autoplot(as.sir(c("S", "I", "R"))), "gg")
}
expect_stdout(print(as.sir(c("S", "I", "R"))))
expect_equal(as.character(as.sir(c(1:3))), c("S", "I", "R"))
expect_equal(as.character(as.sir(c(1:3))), c("S", "I", "R"))
expect_equal(suppressWarnings(as.logical(as.sir("INVALID VALUE"))), NA)
expect_equal(
  summary(as.sir(c("S", "R"))),
  structure(c(
    "Class" = "sir",
    "%R" = "50.0% (n=1)",
    "%SI" = "50.0% (n=1)",
    "- %S" = "50.0% (n=1)",
    "- %I" = " 0.0% (n=0)"
  ), class = c("summaryDefault", "table"))
)
expect_identical(
  as.logical(lapply(example_isolates, is_sir_eligible)),
  as.logical(lapply(example_isolates, is.sir))
)
expect_error(as.sir.mic(as.mic(16)))
expect_error(as.sir.disk(as.disk(16)))
expect_error(get_guideline("this one does not exist"))
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  # 40 sir columns
  expect_equal(
    example_isolates %>%
      mutate_at(vars(PEN:RIF), as.character) %>%
      lapply(is_sir_eligible) %>%
      as.logical() %>%
      sum(),
    40
  )
  expect_equal(sum(is.sir(example_isolates)), 40)

  expect_stdout(print(tibble(ab = as.sir("S"))))
  
  expect_true(example_isolates %>% 
                select(AMC, MEM) %>% 
                mutate(MEM = as.sir(ifelse(AMC == "S", "S", MEM))) %>% 
                pull(MEM) %>% 
                is.sir())
}
if (AMR:::pkg_is_available("skimr", min_version = "2.0.0")) {
  expect_inherits(
    skim(example_isolates),
    "data.frame"
  )
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
    expect_inherits(
      example_isolates %>%
        mutate(
          m = as.mic(2),
          d = as.disk(20)
        ) %>%
        skim(),
      "data.frame"
    )
  }
}

expect_equal(as.sir(c("", "-", NA, "NULL")), c(NA_sir_, NA_sir_, NA_sir_, NA_sir_))

# S. pneumoniae/ampicillin in EUCAST 2020: 0.5-2 ug/ml (R is only > 2)
expect_equal(suppressMessages(
  as.character(
    as.sir(
      x = as.mic(c(0.125, 0.5, 1, 2, 4)),
      mo = "B_STRPT_PNMN",
      ab = "AMP",
      guideline = "EUCAST 2020"
    )
  )),
  c("S", "S", "I", "I", "R")
)
# S. pneumoniae/amoxicillin in CLSI 2019: 2-8 ug/ml (R is 8 and > 8)
expect_equal(suppressMessages(
  as.character(
    as.sir(
      x = as.mic(c(1, 2, 4, 8, 16)),
      mo = "B_STRPT_PNMN",
      ab = "AMX",
      guideline = "CLSI 2019"
    )
  )),
  c("S", "S", "I", "R", "R")
)

expect_true(is.data.frame(sir_interpretation_history(clean = FALSE)))
expect_true(is.data.frame(sir_interpretation_history(clean = TRUE)))
expect_true(is.null(sir_interpretation_history()))

# cutoffs at MIC = 8
expect_equal(
  suppressMessages(as.sir(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2020")),
  as.sir("S")
)
expect_equal(
  suppressMessages(as.sir(as.mic(32), "E. coli", "ampicillin", guideline = "EUCAST 2020")),
  as.sir("R")
)
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  expect_true(suppressWarnings(example_isolates %>%
    mutate(amox_mic = as.mic(2)) %>%
    select(mo, amox_mic) %>%
    as.sir() %>%
    pull(amox_mic) %>%
    is.sir()))
}

expect_equal(
  as.character(
    as.sir(
      x = as.disk(22),
      mo = "B_STRPT_PNMN",
      ab = "ERY",
      guideline = "CLSI"
    )
  ),
  "S"
)
expect_equal(
  as.character(
    as.sir(
      x = as.disk(18),
      mo = "B_STRPT_PNMN",
      ab = "ERY",
      guideline = "CLSI"
    )
  ),
  "I"
)
expect_equal(
  as.character(
    as.sir(
      x = as.disk(10),
      mo = "B_STRPT_PNMN",
      ab = "ERY",
      guideline = "CLSI"
    )
  ),
  "R"
)
if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0")) {
  expect_true(example_isolates %>%
    mutate(amox_disk = as.disk(15)) %>%
    select(mo, amox_disk) %>%
    as.sir(guideline = "CLSI") %>%
    pull(amox_disk) %>%
    is.sir())
  
  # used by group_by() on sir_calc_df(), check some internals to see if grouped calculation without tidyverse works
  groups <- example_isolates %>%
    group_by(mo) %>%
    attributes() %>%
    .$groups
  expect_equal(nrow(groups),
               90)
  expect_equal(class(groups$.rows),
               c("vctrs_list_of", "vctrs_vctr", "list"))
  expect_equal(groups$.rows[[1]],
               c(101, 524, 1368))
  expect_equal(example_isolates[c(101, 524, 1368), "mo", drop = TRUE],
               rep(groups$mo[1], 3))
}
# frequency tables
if (AMR:::pkg_is_available("cleaner")) {
  expect_inherits(cleaner::freq(example_isolates$AMX), "freq")
}


df <- data.frame(
  microorganism = "Escherichia coli",
  AMP = as.mic(8),
  CIP = as.mic(0.256),
  GEN = as.disk(18),
  TOB = as.disk(16),
  ERY = "R", # note about assigning <rsi> class
  CLR = "V"
) # note about cleaning
expect_inherits(
  suppressWarnings(as.sir(df)),
  "data.frame"
)
expect_inherits(
  suppressWarnings(as.sir(data.frame(
    mo = "Escherichia coli",
    amoxi = c("S", "I", "R", "invalid")
  ))$amoxi),
  "sir"
)
expect_warning(as.sir(data.frame(
  mo = "E. coli",
  NIT = c("<= 2", 32)
)))
expect_message(as.sir(data.frame(
  mo = "E. coli",
  NIT = c("<= 2", 32),
  uti = TRUE
)))
expect_message(as.sir(data.frame(
  mo = "E. coli",
  NIT = c("<= 2", 32),
  specimen = c("urine", "blood")
)))

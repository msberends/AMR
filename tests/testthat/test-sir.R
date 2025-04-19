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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

test_that("test-sir.R", {
  skip_on_cran()

  # Existing SIR ------------------------------------------------------------

  # we must only have EUCAST and CLSI, because otherwise the rules in as.sir() will fail
  expect_identical(
    unique(gsub("[^A-Z]", "", AMR::clinical_breakpoints$guideline)),
    c("EUCAST", "CLSI")
  )
  # no missing SDDs
  expect_identical(sum(is.na(AMR::clinical_breakpoints$is_SDD)), 0L)

  expect_true(as.sir("S") < as.sir("I"))
  expect_true(as.sir("I") < as.sir("R"))
  expect_true(is.sir(as.sir("S")))
  x <- example_isolates$AMX
  expect_inherits(x[1], "sir")
  expect_inherits(x[[1]], "sir")
  expect_inherits(c(x[1], x[9]), "sir")
  expect_inherits(unique(x[1], x[9]), "sir")
  pdf(NULL) # prevent Rplots.pdf being created
  expect_silent(barplot(as.sir(c("S", "SDD", "I", "R", "NI"))))
  expect_silent(plot(as.sir(c("S", "SDD", "I", "R", "NI"))))
  if (AMR:::pkg_is_available("ggplot2")) {
    expect_inherits(ggplot2::autoplot(as.sir(c("S", "SDD", "I", "R", "NI"))), "gg")
  }
  expect_output(print(as.sir(c("S", "SDD", "I", "R", "NI"))))
  expect_equal(as.character(as.sir(c(1:3))), c("S", "I", "R"))
  expect_equal(as.character(as.sir(c(1:3))), c("S", "I", "R"))
  expect_equal(suppressWarnings(as.logical(as.sir("INVALID VALUE"))), NA)
  expect_equal(
    summary(as.sir(c("S", "R"))),
    structure(c(
      "Class" = "sir",
      "%S" = "50.0% (n=1)",
      "%SDD" = " 0.0% (n=0)",
      "%I" = " 0.0% (n=0)",
      "%R" = "50.0% (n=1)",
      "%NI" = " 0.0% (n=0)"
    ), class = c("summaryDefault", "table"))
  )
  expect_identical(
    as.logical(lapply(example_isolates, is_sir_eligible)),
    as.logical(lapply(example_isolates, is.sir))
  )
  expect_error(as.sir.mic(as.mic(16)))
  expect_error(as.sir.disk(as.disk(16)))
  expect_error(AMR:::get_guideline("this one does not exist"))

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
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

    expect_output(print(tibble(ab = as.sir("S"))))

    expect_true(example_isolates %>%
      select(AMC, MEM) %>%
      mutate(MEM = as.sir(ifelse(AMC == "S", "S", MEM))) %>%
      pull(MEM) %>%
      is.sir())

    expect_true(example_isolates %>%
      select(AMC, MEM) %>%
      mutate(MEM = if_else(AMC == "S", "S", MEM)) %>%
      pull(MEM) %>%
      is.sir())
  }
  if (AMR:::pkg_is_available("skimr", min_version = "2.0.0", also_load = TRUE)) {
    expect_inherits(
      skim(example_isolates),
      "data.frame"
    )
    if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
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


  # Human -------------------------------------------------------------------

  # allow for guideline length > 1
  expect_equal(
    AMR:::get_guideline(c("CLSI", "CLSI", "CLSI2023", "EUCAST", "EUCAST2020"), AMR::clinical_breakpoints),
    c("CLSI 2024", "CLSI 2024", "CLSI 2023", "EUCAST 2024", "EUCAST 2020")
  )

  # these are used in the script
  expect_true(all(c("B_GRAMN", "B_GRAMP", "B_ANAER-NEG", "B_ANAER-POS", "B_ANAER") %in% AMR::microorganisms$mo))

  mics <- as.mic(2^c(-4:6)) # 0.0625 to 64 in factors of 2
  expect_identical(
    as.character(as.sir(mics,
      mo = "Enterobacterales", ab = "AMC", guideline = "EUCAST 2022",
      uti = FALSE, include_PKPD = FALSE
    )),
    c("S", "S", "S", "S", "S", "S", "S", "S", "R", "R", "R")
  )
  expect_identical(
    as.character(as.sir(mics,
      mo = "Enterobacterales", ab = "AMC", guideline = "EUCAST 2022",
      uti = TRUE, include_PKPD = FALSE
    )),
    c("S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "R")
  )
  expect_identical(
    as.character(as.sir(mics,
      mo = "Escherichia coli", ab = "AMC", guideline = "EUCAST 2022",
      uti = FALSE, include_PKPD = FALSE
    )),
    c("S", "S", "S", "S", "S", "S", "S", "S", "R", "R", "R")
  )

  # test SIR using dplyr's mutate_if(...) and mutate(across(...))
  out1 <- as.sir(as.mic(c(0.256, 0.5, 1, 2)), mo = "Escherichia coli", ab = "ertapenem", guideline = "EUCAST 2023")
  expect_identical(out1, as.sir(c("S", "S", "R", "R")))
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    out2 <- data.frame(
      mo = "Escherichia coli",
      ab = "ertapenem",
      some_mics = as.mic(c(0.256, 0.5, 1, 2))
    ) %>%
      mutate(across(where(is.mic), function(x) as.sir(x, mo = "mo", ab = "ab", guideline = "EUCAST 2023"))) %>%
      pull(some_mics)
    out3 <- data.frame(
      mo = "Escherichia coli",
      ab = "ertapenem",
      some_mics = as.mic(c(0.256, 0.5, 1, 2))
    ) %>%
      mutate_if(is.mic, as.sir, mo = "mo", ab = "ab", guideline = "EUCAST 2023") %>%
      pull(some_mics)

    expect_identical(out1, out2)
    expect_identical(out1, out3)
  }

  # S. pneumoniae/ampicillin in EUCAST 2020: 0.5-2 ug/ml (R is only > 2)
  expect_equal(
    suppressMessages(
      as.character(
        as.sir(
          x = as.mic(c(0.125, 0.5, 1, 2, 4)),
          mo = "B_STRPT_PNMN",
          ab = "AMP",
          guideline = "EUCAST 2020"
        )
      )
    ),
    c("S", "S", "I", "I", "R")
  )
  # S. pneumoniae/amoxicillin in CLSI 2019: 2-8 ug/ml (R is 8 and > 8)
  expect_equal(
    suppressMessages(
      as.character(
        as.sir(
          x = as.mic(c(1, 2, 4, 8, 16)),
          mo = "B_STRPT_PNMN",
          ab = "AMX",
          guideline = "CLSI 2019"
        )
      )
    ),
    c("S", "S", "I", "R", "R")
  )

  expect_true(is.data.frame(sir_interpretation_history(clean = FALSE)))
  expect_true(is.data.frame(sir_interpretation_history(clean = TRUE)))
  expect_true(NROW(sir_interpretation_history()) == 0)

  # cutoffs at MIC = 8
  expect_equal(
    suppressMessages(as.sir(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2020")),
    as.sir("S")
  )
  expect_equal(
    suppressMessages(as.sir(as.mic(32), "E. coli", "ampicillin", guideline = "EUCAST 2020")),
    as.sir("R")
  )
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
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
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
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
    expect_equal(
      nrow(groups),
      90
    )
    expect_equal(
      class(groups$.rows),
      c("vctrs_list_of", "vctrs_vctr", "list")
    )
    expect_equal(
      groups$.rows[[1]],
      c(101, 524, 1368)
    )
    expect_equal(
      example_isolates[c(101, 524, 1368), "mo", drop = TRUE],
      rep(groups$mo[1], 3)
    )
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
      amoxi = c("S", "SDD", "I", "R", "NI", "invalid")
    ))$amoxi),
    "sir"
  )
  # expect_warning(as.sir(data.frame(mo = "E. coli", NIT = c("<= 2", 32))))
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

  # SDD vs I in CLSI 2024
  expect_identical(
    as.sir(as.mic(2^c(-2:4)), mo = "Enterococcus faecium", ab = "Dapto", guideline = "CLSI 2024"),
    as.sir(c("SDD", "SDD", "SDD", "SDD", "SDD", "R", "R"))
  )
  expect_identical(
    as.sir(as.mic(2^c(-2:2)), mo = "Enterococcus faecium", ab = "Cipro
                        ", guideline = "CLSI 2024"),
    as.sir(c("S", "S", "S", "I", "R"))
  )


  # Veterinary --------------------------------------------------------------

  # multiple guidelines
  sir_history <- sir_interpretation_history(clean = TRUE)
  x <- as.sir(as.mic(c(16, 16)), mo = "B_STRPT_CANS", ab = "AMK", host = "dog", guideline = c("CLSI 2024", "CLSI 2014"))
  expect_equal(x, as.sir(c("R", NA)))
  sir_history <- sir_interpretation_history(clean = TRUE)
  expect_equal(sir_history$guideline, c("CLSI 2024", "CLSI 2014"))
  sir_history <- sir_interpretation_history(clean = TRUE)

  mics <- as.mic(2^c(-4:6)) # 0.0625 to 64 in factors of 2
  vet <- data.frame(
    animal = c(rep("cat", 3), rep("dogs", 3), "canine", "equine", "horse", "cattle", "bird"),
    PRA = mics,
    FLR = mics,
    mo = mo_name(rep(c("B_ESCHR_COLI", "B_PSTRL_MLTC", "B_MNNHM_HMLY"), 4)[-1])
  )

  out_vet <- as.sir(vet, host = vet$animal, guideline = "CLSI 2023")
  # give host column name instead of values
  expect_identical(
    out_vet,
    as.sir(vet, host = "animal", guideline = "CLSI 2023")
  )

  # check outcomes
  expect_identical(out_vet$PRA, as.sir(c("S", NA, "S", NA, NA, "R", NA, NA, NA, "R", NA)))
  expect_identical(out_vet$FLR, as.sir(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "R", NA)))

  out_vet <- as.sir(vet, host = "animal", guideline = "EUCAST 2023")
  expect_identical(out_vet$PRA, rep(NA_sir_, 11))
  # expect_identical(out_vet$FLR, as.sir(c("S", "S", NA, "S", "S", NA, "I", "R", NA, "R", "R")))
  expect_identical(out_vet$FLR, as.sir(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "R", NA)))

  sir_history <- sir_interpretation_history()
  expect_identical(
    sort(sir_history$host),
    c(
      "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cats", "cattle",
      "cattle", "cattle", "cattle", "cattle", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs",
      "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "dogs", "horse", "horse", "horse", "horse", "horse", "horse", "horse", "horse",
      "horse", "horse", "poultry", "poultry", "poultry", "poultry", "poultry"
    )
  )

  # ECOFF -------------------------------------------------------------------

  expect_equal(
    suppressMessages(as.sir(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2020", breakpoint_type = "ECOFF")),
    as.sir("S")
  )
  # old method
  expect_warning(as.sir(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2020", ecoff = TRUE))
})

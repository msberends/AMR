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

  # Existing SIR ----------------------------------------------------------

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

  # skimr
  if (AMR:::pkg_is_available("skimr", min_version = "2.0.0", also_load = TRUE)) {
    expect_named(
      skim(example_isolates$PEN),
      c("skim_type", "skim_variable", "n_missing", "complete_rate", "sir.count_S", "sir.count_I", "sir.count_R", "sir.prop_S", "sir.prop_I", "sir.prop_R", "sir.hist")
    )
  }

  expect_equal(as.sir(c("", "-", NA, "NULL")), c(NA_sir_, NA_sir_, NA_sir_, NA_sir_))


  # Human -----------------------------------------------------------------

  # allow for guideline length > 1
  expect_equal(
    AMR:::get_guideline(c("CLSI", "CLSI", "CLSI2023", "EUCAST", "EUCAST2020"), AMR::clinical_breakpoints),
    c("CLSI 2026", "CLSI 2026", "CLSI 2023", "EUCAST 2026", "EUCAST 2020")
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
  ), info = TRUE))
  expect_message(as.sir(data.frame(
    mo = "E. coli",
    NIT = c("<= 2", 32),
    specimen = c("urine", "blood")
  ), info = TRUE))

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


  # Veterinary ------------------------------------------------------------

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

  out_vet <- suppressWarnings(as.sir(vet, host = vet$animal, guideline = "CLSI 2023"))
  # give host column name instead of values
  expect_identical(
    out_vet,
    suppressWarnings(as.sir(vet, host = "animal", guideline = "CLSI 2023"))
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

  # ECOFF -----------------------------------------------------------------

  expect_equal(
    suppressMessages(as.sir(as.mic(c(2, 32)), "E. coli", "ampicillin", guideline = "EUCAST 2020", breakpoint_type = "ECOFF")),
    as.sir(c("WT", "NWT")) # since ECOFF returns WT/NWT at default
  )
  expect_equal(
    suppressMessages(as.sir(as.mic(c(2, 32)), "E. coli", "ampicillin", guideline = "EUCAST 2020", breakpoint_type = "ECOFF", as_wt_nwt = FALSE)),
    as.sir(c("S", "R"))
  )
  # old method
  expect_warning(as.sir(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2020", ecoff = TRUE))


  # Capped MIC handling ---------------------------------------------------

  out1 <- as.sir(as.mic(c("0.125", "<0.125", ">0.125")), mo = "E. coli", ab = "Cipro", guideline = "EUCAST 2025", breakpoint_type = "ECOFF", capped_mic_handling = "none")
  out2 <- as.sir(as.mic(c("0.125", "<0.125", ">0.125")), mo = "E. coli", ab = "Cipro", guideline = "EUCAST 2025", breakpoint_type = "ECOFF", capped_mic_handling = "conservative")
  out3 <- as.sir(as.mic(c("0.125", "<0.125", ">0.125")), mo = "E. coli", ab = "Cipro", guideline = "EUCAST 2025", breakpoint_type = "ECOFF", capped_mic_handling = "standard")
  out4 <- as.sir(as.mic(c("0.125", "<0.125", ">0.125")), mo = "E. coli", ab = "Cipro", guideline = "EUCAST 2025", breakpoint_type = "ECOFF", capped_mic_handling = "lenient")
  expect_equal(out1, as.sir(c("NWT", "NWT", "NWT")))
  expect_equal(out2, as.sir(c("NWT", "NI", "NWT")))
  expect_equal(out3, as.sir(c("NWT", "WT", "NWT")))
  expect_equal(out4, as.sir(c("NWT", "WT", "NWT")))

  # Issue #278: re-running as.sir() on already-<sir> data must preserve columns
  df_already_sir <- data.frame(
    mo  = "B_ESCHR_COLI",
    AMC = as.mic(c("1", "2", "4")),
    GEN = sample(c("S", "I", "R"), 3, replace = TRUE),
    stringsAsFactors = FALSE
  )
  first_pass  <- suppressMessages(as.sir(df_already_sir, col_mo = "mo", info = FALSE))
  second_pass <- suppressMessages(as.sir(first_pass,    col_mo = "mo", info = FALSE))
  expect_equal(ncol(first_pass), ncol(second_pass))
  expect_true(is.sir(second_pass[["AMC"]]))
  expect_true(is.sir(second_pass[["GEN"]]))
  expect_identical(first_pass[["AMC"]], second_pass[["AMC"]])
  expect_identical(first_pass[["GEN"]], second_pass[["GEN"]])

  # Issue #278: metadata columns whose names coincidentally match antibiotic
  # codes (e.g. 'patient' -> OXY, 'ward' -> PRU) must not be processed
  df_meta <- data.frame(
    mo      = "B_ESCHR_COLI",
    patient = paste0("Pt_", 1:20),
    ward    = rep(c("ICU", "Surgery", "Outpatient", "ED"), 5),
    AMC     = as.mic(rep(c("1", "2", "4", "8"), 5)),
    stringsAsFactors = FALSE
  )
  df_meta_sir <- suppressMessages(as.sir(df_meta, col_mo = "mo", info = FALSE))
  expect_true("patient" %in% colnames(df_meta_sir))
  expect_true("ward"    %in% colnames(df_meta_sir))
  expect_false(is.sir(df_meta_sir[["patient"]]))
  expect_false(is.sir(df_meta_sir[["ward"]]))
  expect_true(is.sir(df_meta_sir[["AMC"]]))

  # Parallel computing ----------------------------------------------------
  # Tests must pass even when only 1 core is available; parallel = TRUE then
  # silently falls back to sequential, but results must still be identical.

  set.seed(42)
  n_par <- 200
  df_par <- data.frame(
    mo  = "B_ESCHR_COLI",
    AMC = as.mic(sample(c("0.25", "0.5", "1", "2", "4", "8", "16", "32"), n_par, TRUE)),
    GEN = as.mic(sample(c("0.5", "1", "2", "4", "8", "16", "32", "64"), n_par, TRUE)),
    CIP = as.mic(sample(c("0.001", "0.002", "0.004", "0.008", "0.016", "0.032"), n_par, TRUE)),
    PEN = sample(c("S", "I", "R", NA_character_), n_par, TRUE),
    stringsAsFactors = FALSE
  )

  # clear any existing history before comparing
  sir_interpretation_history(clean = TRUE)
  sir_seq <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE))
  log_seq <- sir_interpretation_history(clean = TRUE)

  sir_par <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE))
  log_par <- sir_interpretation_history(clean = TRUE)

  # 1. parallel = TRUE gives identical SIR results to sequential
  expect_identical(sir_seq[["AMC"]], sir_par[["AMC"]])
  expect_identical(sir_seq[["GEN"]], sir_par[["GEN"]])
  expect_identical(sir_seq[["CIP"]], sir_par[["CIP"]])
  expect_identical(sir_seq[["PEN"]], sir_par[["PEN"]])

  # 2. same number of log rows as sequential
  expect_equal(nrow(log_seq), nrow(log_par))

  # 3. pre-existing log entries must not be duplicated
  #    run sequential once to populate the history, then run parallel and
  #    verify the new parallel run adds exactly as many rows as sequential
  sir_interpretation_history(clean = TRUE)
  suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE)) # populate history
  pre_n <- nrow(sir_interpretation_history())
  suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE))
  post_n <- nrow(sir_interpretation_history())
  expect_equal(post_n - pre_n, nrow(log_seq)) # exactly one run's worth of new rows
  sir_interpretation_history(clean = TRUE)

  # 4. two sequential runs and two parallel runs yield identical results
  sir_par2 <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE))
  expect_identical(sir_par[["AMC"]], sir_par2[["AMC"]])
  expect_identical(sir_par[["GEN"]], sir_par2[["GEN"]])

  # 5. max_cores = 1 gives same results as default sequential
  sir_mc1 <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE, max_cores = 1L))
  expect_identical(sir_seq[["AMC"]], sir_mc1[["AMC"]])
  expect_identical(sir_seq[["GEN"]], sir_mc1[["GEN"]])

  # 6. max_cores = 2 and max_cores = 3 give same results as sequential
  sir_mc2 <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE, max_cores = 2L))
  sir_mc3 <- suppressMessages(as.sir(df_par, col_mo = "mo", info = FALSE, parallel = TRUE, max_cores = 3L))
  expect_identical(sir_seq[["AMC"]], sir_mc2[["AMC"]])
  expect_identical(sir_seq[["GEN"]], sir_mc3[["GEN"]])

  # 7. single-column data frame falls back silently to sequential
  df_single <- df_par[, c("mo", "AMC")]
  sir_single_seq <- suppressMessages(as.sir(df_single, col_mo = "mo", info = FALSE))
  sir_single_par <- suppressMessages(as.sir(df_single, col_mo = "mo", info = FALSE, parallel = TRUE))
  expect_identical(sir_single_seq[["AMC"]], sir_single_par[["AMC"]])

  # 9. row-batch mode (n_cols < n_cores): force row splitting via max_cores and
  #    verify identical output to sequential for a dataset with 2 AB columns so
  #    pieces_per_col = ceiling(max_cores / 2) >= 2 and row batching activates
  df_wide <- data.frame(
    mo  = "B_ESCHR_COLI",
    AMC = as.mic(sample(c("1", "2", "4", "8"), n_par, TRUE)),
    GEN = as.mic(sample(c("1", "2", "4", "8"), n_par, TRUE)),
    stringsAsFactors = FALSE
  )
  sir_wide_seq <- suppressMessages(as.sir(df_wide, col_mo = "mo", info = FALSE))
  sir_wide_par <- suppressMessages(as.sir(df_wide, col_mo = "mo", info = FALSE,
                                          parallel = TRUE, max_cores = 8L))
  expect_identical(sir_wide_seq[["AMC"]], sir_wide_par[["AMC"]])
  expect_identical(sir_wide_seq[["GEN"]], sir_wide_par[["GEN"]])

  # 8. info = TRUE with parallel does not produce per-column worker messages
  #    (messages should only appear in the main process, not duplicated from workers)
  msgs <- capture.output(
    suppressWarnings(as.sir(df_par, col_mo = "mo", info = TRUE, parallel = TRUE)),
    type = "message"
  )
  # each AB column name should appear at most once in all messages combined
  for (ab_nm in c("AMC", "GEN", "CIP", "PEN")) {
    n_mentions <- sum(grepl(ab_nm, msgs, fixed = TRUE))
    expect_lte(n_mentions, 1L)
  }
})

# issue #239 — custom reference_data support
test_that("custom reference_data: non-EUCAST/CLSI guideline produces R", {
  # Build a minimal one-row custom breakpoint table from a plain data.frame.
  # coerce_reference_data_columns() will coerce mo/ab to the right class.
  my_bp <- clinical_breakpoints[clinical_breakpoints$method == "MIC" &
    clinical_breakpoints$type == "human", ][1, ]
  my_bp$guideline    <- "MyLab 2025"
  my_bp$mo           <- "B_ACHRMB_XYLS"  # plain character — coerced to <mo>
  my_bp$ab           <- "MEM"             # plain character — coerced to <ab>
  my_bp$breakpoint_S <- 8
  my_bp$breakpoint_R <- 32

  # guideline omitted: all rows in reference_data are used; R via open interval (>)
  expect_equal(as.character(suppressMessages(
    as.sir(as.mic(64), mo = "B_ACHRMB_XYLS", ab = "MEM", reference_data = my_bp)
  )), "R")
  expect_equal(as.character(suppressMessages(
    as.sir(as.mic(16), mo = "B_ACHRMB_XYLS", ab = "MEM", reference_data = my_bp)
  )), "I")
  # at R breakpoint value must be I (open interval: > not >=)
  expect_equal(as.character(suppressMessages(
    as.sir(as.mic(32), mo = "B_ACHRMB_XYLS", ab = "MEM", reference_data = my_bp)
  )), "I")

  # guideline explicitly set: same result when it matches the data
  expect_equal(as.character(suppressMessages(
    as.sir(as.mic(64), mo = "B_ACHRMB_XYLS", ab = "MEM",
      guideline = "MyLab 2025", reference_data = my_bp)
  )), "R")
})

test_that("custom reference_data: host = NA acts as host-agnostic fallback", {
  my_bp <- clinical_breakpoints[clinical_breakpoints$method == "MIC" &
    clinical_breakpoints$type == "human", ][1, ]
  my_bp$guideline    <- "MyLab 2025"
  my_bp$mo           <- "B_ACHRMB_XYLS"
  my_bp$ab           <- "MEM"
  my_bp$type         <- "animal"
  my_bp$host         <- NA  # logical NA — coerced to character by coerce_reference_data_columns()
  my_bp$breakpoint_S <- 8
  my_bp$breakpoint_R <- 32

  # NA host should match when no species-specific row exists
  result <- suppressMessages(
    as.sir(as.mic(64), mo = "B_ACHRMB_XYLS", ab = "MEM",
      host = "dogs", breakpoint_type = "animal", reference_data = my_bp)
  )
  expect_equal(as.character(result), "R")
})

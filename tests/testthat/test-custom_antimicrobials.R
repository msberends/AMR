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
# Center Groningen in the Netherlands, in collaboration with many      #
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

test_that("test-custom ab.R", {
  skip_on_cran()

  ab_reset_session()

  expect_message(as.ab("testab", info = TRUE))

  suppressMessages(
    add_custom_antimicrobials(
      data.frame(
        ab = "TESTAB",
        name = "Test Antibiotic",
        group = "Test Group"
      )
    )
  )

  expect_identical(as.character(as.ab("testab")), "TESTAB")
  expect_identical(ab_name("testab"), "Test Antibiotic")
  expect_identical(ab_group("testab"), "Test Group")
})

test_that("test-custom ab synonyms", {
  skip_on_cran()

  suppressMessages(clear_custom_antimicrobials())
  suppressMessages(ab_reset_session())
  on.exit(suppressMessages(clear_custom_antimicrobials()), add = TRUE)

  # codes must exist
  expect_error(add_custom_antimicrobial_synonyms("NOTANAB", "something"))

  # ASCII and non-ASCII synonyms (a Korean trade name of piperacillin/tazobactam)
  tazocin_ko <- "\uD0C0\uC870\uC2E0\uC8FC"
  suppressMessages(add_custom_antimicrobial_synonyms("TZP", c("Tazocinum Testname", tazocin_ko)))
  expect_identical(as.character(as.ab("Tazocinum Testname")), "TZP")
  expect_identical(as.character(as.ab(tazocin_ko)), "TZP")
  expect_identical(as.character(as.ab("\uD0C0\uC870 \uC2E0\uC8FC")), "TZP") # white space is ignored
  expect_identical(as.character(as.ab(c(tazocin_ko, "amoxicillin", NA))), c("TZP", "AMX", NA))
  expect_identical(ab_name(tazocin_ko), ab_name("TZP"))

  # another non-ASCII name that was not added must not be matched to TZP
  expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab("\uBA54\uB85C\uD39C")))), "TZP"))

  # data.frame input with column "synonyms", also for custom antimicrobials, and also as a list column
  suppressMessages(add_custom_antimicrobials(data.frame(ab = "TESTSYN", name = "Test Synonym Antibiotic")))
  suppressMessages(add_custom_antimicrobial_synonyms(data.frame(ab = "TESTSYN", synonyms = "Testosyn")))
  expect_identical(as.character(as.ab("testosyn")), "TESTSYN")
  syn_list <- data.frame(ab = "TESTSYN", stringsAsFactors = FALSE)
  syn_list$synonyms <- list(c("Testosyn Two", "Testosyn Three"))
  suppressMessages(add_custom_antimicrobial_synonyms(syn_list))
  expect_identical(as.character(as.ab(c("Testosyn Two", "Testosyn Three"))), c("TESTSYN", "TESTSYN"))

  # one synonym cannot refer to two antimicrobials, nor be the code or name of another antimicrobial
  expect_error(add_custom_antimicrobial_synonyms("MEM", tazocin_ko))
  expect_error(add_custom_antimicrobial_synonyms("TZP", "MEM"))
  expect_error(add_custom_antimicrobial_synonyms("TZP", "Meropenem"))
  # nor take over another antimicrobial's synonym, abbreviation, ATC code, CID or LOINC code
  taken <- "already\\s+identify\\s+another\\s+antimicrobial"
  expect_error(add_custom_antimicrobial_synonyms("MEM", "Zosyn"), taken)
  expect_error(add_custom_antimicrobial_synonyms("MEM", "J01CR05"), taken)
  expect_identical(as.character(as.ab(c("Zosyn", "J01CR05"))), c("TZP", "TZP"))
  tzp_row <- AMR::antimicrobials[AMR::antimicrobials$ab == "TZP", , drop = FALSE]
  tzp_ids <- c(
    unlist(tzp_row$synonyms)[1], unlist(tzp_row$abbreviations)[1], unlist(tzp_row$atc)[1],
    as.character(tzp_row$cid[1]), unlist(tzp_row$loinc)[1]
  )
  for (id in tzp_ids[!is.na(tzp_ids) & tzp_ids != ""]) {
    expect_error(add_custom_antimicrobial_synonyms("MEM", id), taken, info = id)
  }
  # also not in the spelling variants that as.ab() already resolves to another antimicrobial
  expect_error(add_custom_antimicrobial_synonyms("MEM", "Piperacillin-tazobactam"), taken)
  expect_error(add_custom_antimicrobial_synonyms("MEM", "Z\u00F3syn"), taken)
  # nor with other punctuation (including an invisible soft hyphen)
  for (id in c("J01-CR05", "J01.CR05", "Zo-syn", "Zo\u00ADsyn")) {
    expect_error(add_custom_antimicrobial_synonyms("MEM", id), taken, info = id)
  }
  # but an antimicrobial may get its own identifiers again
  expect_message(add_custom_antimicrobial_synonyms("TZP", c("Zosyn", "J01CR05")))
  # synonyms consisting of white space only are skipped
  expect_message(add_custom_antimicrobial_synonyms("TZP", c("\u00A0", "\u3000")), "No synonyms")

  # all Unicode white space is ignored (no-break space, ideographic space)
  expect_identical(as.character(as.ab("\uD0C0\uC870\u3000\uC2E0\uC8FC")), "TZP")
  expect_identical(as.character(as.ab("\uD0C0\uC870\u00A0\uC2E0\uC8FC")), "TZP")
  expect_identical(as.character(as.ab("Tazocinum\u00A0Testname")), "TZP")

  # a trailing part in parentheses or a strength is removed when there is no exact match
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " 4.5g"))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, "4.5g"))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " 4.5g "))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " 4.5\u00A0g"))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " 4/0.5 g"))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " 500 mg/mL"))), "TZP")
  expect_identical(as.character(as.ab(paste0(tazocin_ko, " (piperacillin, tazobactam) 4.5g"))), "TZP")
  expect_identical(as.character(as.ab("Tazocinum Testname 2.25 g")), "TZP")
  # also nested parentheses, brackets, trailing punctuation, unit symbols and full-width characters
  tazocin_variants <- paste0(tazocin_ko, c(
    "4.5g(pip(4:1))", " [4.5g]", "-4.5g", " 4.5g.", " 4.5\u338E", " \uFF14\uFF0E\uFF15\uFF47"
  ))
  expect_identical(as.character(as.ab(tazocin_variants)), rep("TZP", length(tazocin_variants)))
  # but nothing is matched on a strength alone
  expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab("4.5g")))), "TZP"))
  # and no words are removed, in any language, such as dosage forms or pack quantities
  for (word in c(" Inj", " IV", " 1 vial", "\uC8FC\uC0AC", " 1\uBCD1")) {
    expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab(paste0("Tazocinum Testname", word))))), "TZP"), info = word)
  }

  # one trailing part is removed at a time, so the longest registered synonym wins
  suppressMessages(add_custom_antimicrobial_synonyms("CXA", "Cefubrand Testname"))
  suppressMessages(add_custom_antimicrobial_synonyms("CXM", "Cefubrand Testname (IV)"))
  expect_identical(
    as.character(as.ab(c("Cefubrand Testname (IV) 750mg", "Cefubrand Testname (IV) 750 mg.", "Cefubrand Testname 250mg"))),
    c("CXM", "CXM", "CXA")
  )

  # synonyms are listed by ab_synonyms(), but take no part in fuzzy matching
  expect_true("Tazocinum Testname" %in% ab_synonyms("TZP"))
  tzp <- which(AMR:::AMR_env$AB_lookup$ab == "TZP")
  expect_false(AMR:::generalise_antibiotic_name("Tazocinum Testname") %in% AMR:::AMR_env$AB_lookup$generalised_all[[tzp]])
  expect_false(AMR:::generalise_antibiotic_name("Tazocinum Testname") %in% AMR:::AMR_env$AB_lookup$generalised_synonyms[[tzp]])

  # input that is not valid UTF-8, such as a CP949 export in a UTF-8 locale, must not give an error
  tazocin_cp949 <- tryCatch(iconv(tazocin_ko, from = "UTF-8", to = "CP949"), error = function(e) NA_character_)
  if (!is.na(tazocin_cp949) && isTRUE(l10n_info()[["UTF-8"]])) {
    Encoding(tazocin_cp949) <- "unknown"
    expect_true(is.na(AMR:::custom_ab_synonym_key(tazocin_cp949)))
    expect_identical(as.character(suppressWarnings(as.ab(c(tazocin_cp949, tazocin_ko)))), c(NA, "TZP"))
    expect_warning(as.ab(tazocin_cp949))
    expect_error(suppressMessages(add_custom_antimicrobial_synonyms("TZP", tazocin_cp949)), "UTF-8")
  }

  # clearing removes the synonyms as well
  suppressMessages(clear_custom_antimicrobials())
  expect_false(identical(as.character(suppressWarnings(suppressMessages(as.ab(tazocin_ko)))), "TZP"))
})

test_that("test-custom ab synonyms via add_custom_antimicrobials()", {
  skip_on_cran()

  suppressMessages(clear_custom_antimicrobials())
  suppressMessages(ab_reset_session())
  on.exit(suppressMessages(clear_custom_antimicrobials()), add = TRUE)
  tazocin_ko <- "\uD0C0\uC870\uC2E0\uC8FC"

  # rows for existing codes with only "ab" and "synonyms" filled are passed on as synonyms;
  # one synonym per row, codes may repeat, empty cells are skipped (this is what the AMR_custom_ab option loads at start-up)
  syn <- data.frame(
    ab = c("TZP", "TZP", "MEM", "MEM", NA),
    synonyms = c("Tazocinum Testname", tazocin_ko, "Meropenemum Testname", NA, NA),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".rds")
  saveRDS(syn, tmp)
  suppressMessages(add_custom_antimicrobials(readRDS(tmp)))
  unlink(tmp)
  expect_identical(as.character(as.ab(c("Tazocinum Testname", tazocin_ko, "Meropenemum Testname"))), c("TZP", "TZP", "MEM"))

  # existing codes with any other column filled still give the current error
  expect_error(add_custom_antimicrobials(data.frame(ab = "TZP", name = "Something", synonyms = "Something")))
  expect_error(add_custom_antimicrobials(data.frame(ab = "TZP", name = "Something")))

  # plain records still work as before, also without a "synonyms" column
  suppressMessages(add_custom_antimicrobials(data.frame(ab = "TESTREC", name = "Test Record Antibiotic")))
  expect_identical(ab_name("TESTREC"), "Test Record Antibiotic")

  # new antimicrobials and synonyms for them can be added in one go
  suppressMessages(add_custom_antimicrobials(data.frame(
    ab = c("TESTMIX", "TESTMIX"),
    name = c("Test Mix Antibiotic", NA),
    synonyms = c(NA, "Testomix"),
    stringsAsFactors = FALSE
  )))
  expect_identical(as.character(as.ab("Testomix")), "TESTMIX")
  expect_identical(ab_name("Testomix"), "Test Mix Antibiotic")

  # an error in the synonyms leaves the session unchanged, also for new antimicrobials in the same call
  n_before <- NROW(AMR:::AMR_env$AB_lookup)
  expect_error(add_custom_antimicrobials(data.frame(ab = c("NOTANAB2", "TESTNEW"), name = c(NA, "Test New"), synonyms = c("Notanab", NA), stringsAsFactors = FALSE)))
  expect_error(add_custom_antimicrobials(data.frame(ab = c("TESTNEW", "MEM"), name = c("Test New", NA), synonyms = c(NA, "Tazocinum Testname"), stringsAsFactors = FALSE)))
  expect_identical(NROW(AMR:::AMR_env$AB_lookup), n_before)
  expect_false("TESTNEW" %in% AMR:::AMR_env$AB_lookup$ab)

  # a new antimicrobial and a synonym of another one in the same call must not share a name
  expect_error(add_custom_antimicrobials(data.frame(ab = c("TESTSAME", "TZP"), name = c("Foobarcillin", NA), synonyms = c(NA, "Foobarcillin"), stringsAsFactors = FALSE)), "already\\s+identify\\s+another\\s+antimicrobial")
  expect_false("TESTSAME" %in% AMR:::AMR_env$AB_lookup$ab)
  # and a new antimicrobial cannot take a name that is already a user-added synonym
  suppressMessages(add_custom_antimicrobial_synonyms("TZP", "Foobarcillin"))
  expect_error(add_custom_antimicrobials(data.frame(ab = "TESTLATER", name = "Foobarcillin")), "already\\s+in\\s+use\\s+as\\s+user-added\\s+synonyms")
  expect_identical(as.character(as.ab("Foobarcillin")), "TZP")
  # not with other punctuation either, so that the order of adding does not matter
  suppressMessages(add_custom_antimicrobial_synonyms("TZP", "Foo-quxcillin"))
  expect_error(add_custom_antimicrobials(data.frame(ab = "TESTORDER", name = "Foo quxcillin")), "already\\s+in\\s+use\\s+as\\s+user-added\\s+synonyms")

  # a new record without a name does not block the synonyms in the same call
  suppressMessages(add_custom_antimicrobials(data.frame(
    ab = c("TESTNONAME", "MEM"), name = c(NA, NA), group = c("Test", NA), synonyms = c(NA, "Meropenemum Somebrand"),
    stringsAsFactors = FALSE
  )))
  expect_identical(as.character(as.ab("Meropenemum Somebrand")), "MEM")
})

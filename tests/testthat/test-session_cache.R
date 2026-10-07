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

test_that("test-session_cache.R", {
  skip_on_cran()

  # as.ab() ----
  ab_reset_session()
  # a `fast_mode` result must never be reused by a regular coercion
  expect_true(is.na(suppressWarnings(as.ab("amxcilin", fast_mode = TRUE, info = FALSE))))
  expect_identical(as.character(suppressWarnings(as.ab("amxcilin", info = FALSE))), "AMX")
  expect_true(all(c("amxcilin") %in% AMR:::AMR_env$ab_previously_coerced$input))
  expect_identical(colnames(AMR:::AMR_env$ab_previously_coerced), c("key", "input", "value"))
  # second run comes from the cache and gives the same result
  expect_identical(as.character(as.ab("amxcilin", info = FALSE)), "AMX")
  # ATC codes and synonyms are looked up in the session index
  expect_identical(as.character(as.ab(c("J01CA04", "J01MA02", "amoxil", "cipro"), info = FALSE)), c("AMX", "CIP", "AMX", "CIP"))
  # resetting does not touch microorganism uncertainties
  suppressWarnings(as.mo("Klebsiela pneumonia", info = FALSE))
  n_uncertain <- NROW(mo_uncertainties())
  expect_message(ab_reset_session())
  expect_identical(NROW(mo_uncertainties()), n_uncertain)
  expect_identical(NROW(AMR:::AMR_env$ab_previously_coerced), 0L)

  # as.mo() ----
  mo_reset_session()
  x <- c("Klebsiela pneumonia", "Stafylococcus aureus")
  expect_identical(as.character(suppressWarnings(as.mo(x, info = FALSE))), c("B_KLBSL_PNMN", "B_STPHY_AURS"))
  expect_true(all(x %in% AMR:::AMR_env$mo_previously_coerced$input))
  expect_identical(NROW(mo_uncertainties()), 2L)
  # second run comes from the cache, and the uncertainties are still complete
  expect_identical(as.character(suppressWarnings(as.mo(rev(x), info = FALSE))), c("B_STPHY_AURS", "B_KLBSL_PNMN"))
  expect_identical(sort(mo_uncertainties()$original_input), sort(x))
  # another minimum_matching_score is another cache key
  n_cached <- NROW(AMR:::AMR_env$mo_previously_coerced)
  suppressWarnings(as.mo(x, minimum_matching_score = 0.1, info = FALSE))
  expect_identical(NROW(AMR:::AMR_env$mo_previously_coerced), n_cached + 2L)
  # cached results are keyed on the original input, also when it is cleaned before fuzzy matching
  mo_reset_session()
  suppressWarnings(as.mo("Klebsiela pneumonia group", info = FALSE))
  expect_true("Klebsiela pneumonia group" %in% AMR:::AMR_env$mo_previously_coerced$input)
  # input that could not be coerced is also remembered
  expect_identical(as.character(suppressWarnings(as.mo("xq", info = FALSE))), "UNKNOWN")
  expect_warning(as.mo("xq", info = FALSE))
  expect_identical(mo_failures(), "xq")

  # valid codes of outdated names are replaced, also on the fast path
  syn <- AMR:::AMR_env$MO_lookup$mo[AMR:::AMR_env$MO_lookup$status == "synonym"][1:5]
  expect_false(any(as.mo(syn, keep_synonyms = FALSE, info = FALSE) %in% syn, na.rm = TRUE))
  expect_identical(as.character(as.mo(syn, keep_synonyms = TRUE, info = FALSE)), as.character(syn))

  # SNOMED codes
  expect_identical(as.character(as.mo(c(112283007, 3092008), info = FALSE)), c("B_ESCHR_COLI", "B_STPHY_AURS"))

  # synonym_mo_to_accepted_mo() ----
  all_mo <- as.character(AMR:::AMR_env$MO_lookup$mo)
  for (fill in c(TRUE, FALSE)) {
    AMR:::reset_mo_cache()
    expect_identical(
      AMR:::synonym_mo_to_accepted_mo(c(all_mo, NA, "invalid"), fill_in_accepted = fill),
      AMR:::synonym_mo_to_accepted_mo_uncached(c(all_mo, NA, "invalid"), fill_in_accepted = fill, dataset = AMR:::AMR_env$MO_lookup)
    )
    # now from the memoised mapping
    expect_identical(
      AMR:::synonym_mo_to_accepted_mo(rev(all_mo), fill_in_accepted = fill),
      AMR:::synonym_mo_to_accepted_mo_uncached(rev(all_mo), fill_in_accepted = fill, dataset = AMR:::AMR_env$MO_lookup)
    )
  }

  # mo_*() functions ----
  x <- as.mo(c("B_ESCHR_COLI", "B_STPHY_AURS", NA, "B_ESCHR_COLI", "UNKNOWN"))
  expect_identical(mo_name(x, language = NULL), c("Escherichia coli", "Staphylococcus aureus", NA, "Escherichia coli", "(unknown name)"))
  expect_identical(mo_genus(rep(x, 3), language = NULL), rep(mo_genus(x, language = NULL), 3))
  expect_identical(mo_snomed(x[c(1, 1, 3)]), list(mo_snomed(x[1])[[1]], mo_snomed(x[1])[[1]], NULL))
  expect_identical(mo_property(x, "prevalence"), as.double(AMR:::AMR_env$MO_lookup$prevalence[match(x, AMR:::AMR_env$MO_lookup$mo)]))
  # input that already contains the property is returned as is
  expect_identical(mo_genus(c("Escherichia", "Escherichia", NA)), c("Escherichia", "Escherichia", NA))
  # character input with duplicates
  expect_identical(mo_name(c("E. coli", "S. aureus", "E. coli"), language = NULL), c("Escherichia coli", "Staphylococcus aureus", "Escherichia coli"))
  # factors
  expect_identical(mo_gramstain(factor(c("E. coli", "S. aureus", "E. coli")), language = NULL), c("Gram-negative", "Gram-positive", "Gram-negative"))
  # synonyms are replaced unless keep_synonyms = TRUE
  expect_false(any(mo_name(syn, keep_synonyms = FALSE, language = NULL) == mo_name(syn, keep_synonyms = TRUE, language = NULL)))

  # invalidation ----
  mo_reset_session()
  suppressWarnings(as.mo("Klebsiela pneumonia", info = FALSE))
  invisible(AMR:::get_mo_index())
  invisible(AMR:::synonym_mo_to_accepted_mo("B_ESCHR_COLI"))
  suppressMessages(add_custom_microorganisms(data.frame(genus = "Testbacterium", species = "cachei")))
  expect_identical(NROW(AMR:::AMR_env$mo_previously_coerced), 0L)
  expect_null(AMR:::AMR_env$MO_index)
  expect_null(AMR:::AMR_env$mo_accepted)
  expect_identical(mo_name("Testbacterium cachei", language = NULL), "Testbacterium cachei")
  suppressMessages(clear_custom_microorganisms())
  expect_identical(NROW(AMR:::AMR_env$mo_previously_coerced), 0L)

  ab_reset_session()
  expect_true(is.na(suppressWarnings(as.ab("Qwxzqwxz", info = FALSE))))
  invisible(AMR:::get_ab_index())
  suppressMessages(add_custom_antimicrobials(data.frame(ab = "TESTCD", name = "Qwxzqwxz", group = "Test group")))
  expect_identical(NROW(AMR:::AMR_env$ab_previously_coerced), 0L)
  expect_null(AMR:::AMR_env$AB_index)
  expect_identical(as.character(as.ab("Qwxzqwxz", info = FALSE)), "TESTCD")
  suppressMessages(clear_custom_antimicrobials())
  expect_true(is.na(suppressWarnings(as.ab("Qwxzqwxz", info = FALSE))))
})

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

# Tests of the `microorganisms` data set and the stability of MO codes.
# The rules themselves are in helper-microorganisms.R, which is also used by the taxonomy build.
#
# KNOWN DEFECTS: the data set on the main branch of October 2026 contains defects that will be solved by the
# next taxonomy rebuild (e.g. Trichophyton as a bacterium, issue #309). These are listed in
# data-raw/microorganisms_files/microorganisms_known_defects.csv, so that these tests fail on every NEW defect,
# while the current ones are visible. This list may only become shorter, and must be empty after the rebuild (then
# delete the file). As data-raw/ does not ship with the package, the tests that use this list only run in the
# source repository.

known_defects <- function() {
  read_mo_known_defects()
}
new_defects <- function(rule, records) {
  setdiff(records, known_defects()$record[known_defects()$rule == rule])
}

# Sentinel organisms: clinically important taxa with properties that must not change unnoticed.
# NA means: not checked.
sentinels <- data.frame(
  fullname = c(
    "Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa", "Acinetobacter baumannii",
    "Haemophilus influenzae", "Neisseria gonorrhoeae", "Staphylococcus aureus", "Streptococcus pneumoniae",
    "Enterococcus faecium", "Clostridioides difficile", "Mycobacterium tuberculosis", "Listeria monocytogenes",
    "Tropheryma whipplei", "Candida albicans", "Candidozyma auris", "Nakaseomyces glabratus",
    "Aspergillus fumigatus", "Cryptococcus neoformans", "Pneumocystis jirovecii", "Trichophyton rubrum",
    "Mucor", "Necator americanus", "Ascaris lumbricoides", "Echinococcus granulosus",
    "Capillaria", "Plasmodium falciparum", "Leishmania donovani", "Trypanosoma cruzi",
    "Giardia", "Toxoplasma gondii"
  ),
  domain = c(
    rep("Bacteria", 13), rep("Fungi", 8), rep("Animalia", 4),
    "Protozoa", "Protozoa", "Protozoa", "Protozoa", NA
  ),
  gramstain = c(
    rep("Gram-negative", 6), rep("Gram-positive", 7), rep(NA, 17)
  ),
  family = c(
    "Enterobacteriaceae", "Enterobacteriaceae", "Pseudomonadaceae", "Moraxellaceae", "Pasteurellaceae",
    "Neisseriaceae", "Staphylococcaceae", "Streptococcaceae", "Enterococcaceae", "Peptostreptococcaceae",
    "Mycobacteriaceae", "Listeriaceae", NA, NA, NA, NA,
    "Aspergillaceae", NA, NA, "Arthrodermataceae", "Mucoraceae", NA, NA, NA,
    NA, NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)
sentinel_defects <- function() {
  i <- match(sentinels$fullname, microorganisms$fullname)
  ok <- !is.na(i)
  # (paste() of an empty vector would return one element, hence the use of sprintf())
  out <- sprintf("%s is missing", sentinels$fullname[!ok])
  check <- function(property, actual) {
    expected <- sentinels[[property]]
    wrong <- which(ok & !is.na(expected) & actual != expected)
    sprintf("%s %s is %s instead of %s", sentinels$fullname[wrong], rep(property, length(wrong)), actual[wrong], expected[wrong])
  }
  out <- c(out, check("domain", microorganisms$domain[i]))
  out <- c(out, check("family", microorganisms$family[i]))
  gramstain <- rep(NA_character_, length(i))
  gramstain[ok] <- suppressWarnings(mo_gramstain(microorganisms$mo[i[ok]], language = NULL))
  out <- c(out, check("gramstain", gramstain))
  # the code prefix must always match the domain
  prefix_ok <- startsWith(as.character(microorganisms$mo[i]), paste0(MO_DOMAIN_PREFIX[sentinels$domain], "_"))
  c(out, sprintf("%s has a code with the wrong prefix", sentinels$fullname[which(ok & !is.na(sentinels$domain) & !prefix_ok)]))
}

test_that("microorganisms: integrity rules", {
  skip_on_cran()
  skip_if_no_mo_repository()
  issues <- mo_integrity_issues(microorganisms, registry = read_mo_registry(), renames = read_mo_renames(), retirements = read_mo_retirements())
  for (rule in names(issues)) {
    defects <- new_defects(rule, issues[[rule]])
    expect_true(length(defects) == 0,
      info = paste0("Rule '", rule, "' is broken by: ", paste(utils::head(defects, 25), collapse = ", "))
    )
  }
})

test_that("microorganisms: known defects are no longer than needed", {
  skip_on_cran()
  skip_if_no_mo_repository()
  # a known defect that is solved must be removed from the list, so that it cannot come back unnoticed
  kd <- known_defects()
  skip_if(nrow(kd) == 0)
  issues <- mo_integrity_issues(microorganisms, registry = read_mo_registry(), renames = read_mo_renames(), retirements = read_mo_retirements())
  current <- c(
    unlist(lapply(names(issues), function(rule) paste(rule, issues[[rule]], sep = "|"))),
    paste("sentinel", sentinel_defects(), sep = "|")
  )
  solved <- setdiff(paste(kd$rule, kd$record, sep = "|"), current)
  expect_true(length(solved) == 0,
    info = paste0("Solved, remove from microorganisms_known_defects.csv: ", paste(utils::head(solved, 25), collapse = ", "))
  )
})

test_that("microorganisms: sentinel organisms", {
  skip_on_cran()
  skip_if_no_mo_repository()
  defects <- new_defects("sentinel", sentinel_defects())
  expect_true(length(defects) == 0, info = paste(defects, collapse = "; "))
})

test_that("MO codes: every code of every release since v2.0.0 still means the same taxon", {
  skip_on_cran()
  registry <- read_mo_registry()
  skip_if(is.null(registry), "MO code registry not available (only in the source repository)")
  renames <- read_mo_renames()

  # the registry itself
  expect_false(anyDuplicated(registry$mo) > 0)
  expect_true(all(numeric_version(sub("^v", "", registry$first_release)) >= "2.0.0"))
  expect_true(all(numeric_version(sub("^v", "", registry$first_release)) <= numeric_version(sub("^v", "", registry$last_release))))

  # every registered code must lead to its own taxon (or to NA if the taxon is not in the data anymore)
  expected_name <- registry$fullname
  renamed <- match(registry$mo, renames$mo)
  expected_name[!is.na(renamed)] <- renames$new_name[renamed[!is.na(renamed)]]
  # deliberate in as.mo() since v2.0.0: the kingdom/domain Fungi returns the code for an unknown fungus
  expected_name[registry$mo == "F_[KNG]_FUNGI"] <- "(unknown fungus)"
  # codes retired on purpose (their taxon was another organism) must lead to NA
  retirements <- read_mo_retirements()
  retired_on_purpose <- registry$mo %in% retirements$mo
  expect_true(all(is.na(suppressWarnings(as.mo(registry$mo[retired_on_purpose], info = FALSE)))))
  registry <- registry[!retired_on_purpose, , drop = FALSE]
  expected_name <- expected_name[!retired_on_purpose]
  result <- suppressWarnings(as.mo(registry$mo, keep_synonyms = TRUE, info = FALSE))
  result_name <- microorganisms$fullname[match(result, microorganisms$mo)]
  wrong <- !is.na(result) & mo_name_without_suffix(result_name) != mo_name_without_suffix(expected_name)
  # (paste() of empty vectors would return one element, hence the use of sprintf())
  defects <- new_defects("registered_code_with_other_taxon", sprintf("%s (%s)", registry$mo[wrong], result_name[wrong]))
  expect_true(length(defects) == 0,
    info = paste0("Registered codes leading to another taxon: ", paste(utils::head(defects, 25), collapse = ", "))
  )

  # the internal lookup table must be up to date with the registry
  full_registry <- read_mo_registry()
  expect_setequal(AMR:::MO_RETIRED_CODES$old_mo, full_registry$mo[!full_registry$mo %in% microorganisms$mo])
})

test_that("MO codes: codes of earlier releases are translated, never guessed", {
  skip_on_cran()
  retired <- AMR:::MO_RETIRED_CODES
  expect_false(any(retired$old_mo %in% microorganisms$mo))
  expect_true(all(retired$mo[!is.na(retired$mo)] %in% microorganisms$mo))
  # every translated code leads to a record with the registered name
  translated <- retired[!is.na(retired$mo), , drop = FALSE]
  expect_identical(
    as.character(suppressWarnings(as.mo(translated$old_mo, keep_synonyms = TRUE, info = FALSE))),
    translated$mo
  )
  expect_identical(
    mo_name_without_suffix(microorganisms$fullname[match(translated$mo, microorganisms$mo)]),
    mo_name_without_suffix(translated$fullname)
  )

  # warnings are wrapped to the console width, so the patterns allow a line break between words
  # regression: in October 2026, the code of Trichophyton (F_TRCHP, v3.0.1) was guessed as Talaromyces; since the
  # taxonomy update of October 2026, F_TRCHP_RBRM is the current code of Trichophyton rubrum again
  x <- suppressWarnings(as.mo("F_TRCHP_RBRM", info = FALSE))
  expect_identical(mo_name(x, language = NULL), "Trichophyton rubrum")
  expect_identical(suppressWarnings(mo_genus("F_TRCHP", language = NULL)), "Trichophyton")

  # codes from before v2.0.0 and unknown codes are never guessed (decision by Matthijs S. Berends, 3 October 2026)
  expect_warning(x <- as.mo("B_ESCH_COL", info = FALSE), "not\\s+supported")
  expect_true(is.na(x))
  expect_warning(x <- as.mo(c("B_ESCHR_COL", "B_ESCHR_COLI"), info = FALSE), "unknown\\s+MO\\s+code")
  expect_identical(as.character(x), c(NA, "B_ESCHR_COLI"))

  # regression: until October 2026, a MycoBank target that was not in the data set replaced a valid GBIF target by NA
  expect_identical(mo_current("Chaetomium abuense"), "Amesia nigricolor")
  # and a synonym without any known current name is not 'renamed' to itself
  no_target <- mo_synonyms_without_current_name(microorganisms)
  expect_true(all(is.na(AMR:::synonym_mo_to_accepted_mo(no_target$mo[is.na(no_target$chain_ends_at)], fill_in_accepted = FALSE))))

  # current codes are returned as they are, without warnings
  expect_silent(x <- as.mo(c("B_ESCHR_COLI", "B_STPHY_AURS", NA), info = FALSE))
  expect_identical(as.character(x), c("B_ESCHR_COLI", "B_STPHY_AURS", NA))
})

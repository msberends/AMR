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

test_that("data works", {
  # IDs should always be unique
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))
  expect_identical(class(microorganisms$mo), c("mo", "character"))
  expect_identical(nrow(antimicrobials), length(unique(AMR::antimicrobials$ab)))
  expect_true(all(is.na(AMR::antimicrobials$atc[duplicated(AMR::antimicrobials$atc)])))
  expect_identical(class(AMR::antimicrobials$ab), c("ab", "character"))


  # check cross table reference
  expect_true(all(microorganisms.codes$mo %in% microorganisms$mo))
  expect_true(all(example_isolates$mo %in% microorganisms$mo))
  expect_true(all(microorganisms.groups$mo %in% microorganisms$mo))
  expect_true(all(microorganisms.groups$mo_group %in% microorganisms$mo))
  expect_true(all(clinical_breakpoints$mo %in% microorganisms$mo))
  expect_true(all(clinical_breakpoints$ab %in% AMR::antimicrobials$ab))
  expect_true(all(intrinsic_resistant$mo %in% microorganisms$mo))
  expect_true(all(intrinsic_resistant$ab %in% AMR::antimicrobials$ab))
  expect_false(any(is.na(microorganisms.codes$code)))
  expect_false(any(is.na(microorganisms.codes$mo)))
  expect_true(all(dosage$ab %in% AMR::antimicrobials$ab))
  expect_true(all(dosage$name %in% AMR::antimicrobials$name))
  # check valid disks/MICs
  expect_false(any(is.na(as.mic(clinical_breakpoints[which(clinical_breakpoints$method == "MIC" & clinical_breakpoints$ref_tbl != "ECOFF"), "breakpoint_S", drop = TRUE]))))
  expect_false(any(is.na(as.mic(clinical_breakpoints[which(clinical_breakpoints$method == "MIC" & clinical_breakpoints$ref_tbl != "ECOFF"), "breakpoint_R", drop = TRUE]))))
  expect_false(any(is.na(as.disk(clinical_breakpoints[which(clinical_breakpoints$method == "DISK" & clinical_breakpoints$ref_tbl != "ECOFF"), "breakpoint_S", drop = TRUE]))))
  expect_false(any(is.na(as.disk(clinical_breakpoints[which(clinical_breakpoints$method == "DISK" & clinical_breakpoints$ref_tbl != "ECOFF"), "breakpoint_R", drop = TRUE]))))

  # antibiotic names must always be coercible to their original AB code
  expect_identical(as.ab(AMR::antimicrobials$name), AMR::antimicrobials$ab)

  if (AMR:::pkg_is_available("tibble")) {
    # there should be no diacritics (i.e. non ASCII) characters in the datasets (CRAN policy)
    datasets <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item", drop = TRUE]
    for (i in seq_len(length(datasets))) {
      dataset <- get(datasets[i], envir = asNamespace("AMR"))
      expect_identical(AMR:::dataset_UTF8_to_ASCII(dataset), dataset, info = datasets[i])
    }
  }

  df <- AMR:::AMR_env$MO_lookup
  expect_true(all(c(
    "mo", "fullname", "status", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "subspecies", "rank", "ref", "source",
    "lpsn", "lpsn_parent", "lpsn_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to", "prevalence",
    "snomed", "kingdom_index", "fullname_lower", "full_first", "species_first"
  ) %in% colnames(df)))

  expect_inherits(AMR:::MO_CONS, "mo")

  uncategorised <- subset(
    microorganisms,
    genus == "Staphylococcus" &
      !species %in% c("", "aureus") &
      !mo %in% c(AMR:::MO_CONS, AMR:::MO_COPS)
  )
  expect_true(NROW(uncategorised) == 0,
    info = ifelse(NROW(uncategorised) == 0,
      "All staphylococcal species categorised as CoNS/CoPS.",
      paste0(
        "Staphylococcal species not categorised as CoNS/CoPS: S. ",
        uncategorised$species, " (", uncategorised$mo, ")",
        collapse = "\n"
      )
    )
  )

  # THIS WILL CHECK NON-ASCII STRINGS IN ALL FILES:

  # check_non_ascii <- function() {
  #   purrr::map_df(
  #     .id = "file",
  #     # list common text files
  #     .x = fs::dir_ls(
  #       recurse = TRUE,
  #       type = "file",
  #       # ignore images, compressed
  #       regexp = "\\.(png|ico|rda|ai|tar.gz|zip|xlsx|csv|pdf|psd)$",
  #       invert = TRUE
  #     ),
  #     .f = function(path) {
  #       x <- readLines(path, warn = FALSE)
  #       # from tools::showNonASCII()
  #       asc <- iconv(x, "latin1", "ASCII")
  #       ind <- is.na(asc) | asc != x
  #       # make data frame
  #       if (any(ind)) {
  #         tibble::tibble(
  #           row = which(ind),
  #           line = iconv(x[ind], "latin1", "ASCII", sub = "byte")
  #         )
  #       } else {
  #         tibble::tibble()
  #       }
  #     }
  #   )
  # }
  # x <- check_non_ascii() %>%
  #   filter(file %unlike% "^(data-raw|docs|git_)")
})

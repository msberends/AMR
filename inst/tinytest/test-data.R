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
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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

# IDs should always be unique
expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))
expect_identical(class(microorganisms$mo), c("mo", "character"))
expect_identical(nrow(antibiotics), length(unique(antibiotics$ab)))
expect_true(all(is.na(antibiotics$atc[duplicated(antibiotics$atc)])))
expect_identical(class(antibiotics$ab), c("ab", "character"))


# check cross table reference
expect_true(all(microorganisms.codes$mo %in% microorganisms$mo))
expect_true(all(example_isolates$mo %in% microorganisms$mo))
expect_true(all(rsi_translation$mo %in% microorganisms$mo))
expect_true(all(rsi_translation$ab %in% antibiotics$ab))
expect_true(all(intrinsic_resistant$mo %in% microorganisms$mo))
expect_true(all(intrinsic_resistant$ab %in% antibiotics$ab))
expect_false(any(is.na(microorganisms.codes$code)))
expect_false(any(is.na(microorganisms.codes$mo)))
expect_true(all(dosage$ab %in% antibiotics$ab))
expect_true(all(dosage$name %in% antibiotics$name))
# check valid disks/MICs
expect_false(any(is.na(as.mic(rsi_translation[which(rsi_translation$method == "MIC"), "breakpoint_S", drop = TRUE]))))
expect_false(any(is.na(as.mic(rsi_translation[which(rsi_translation$method == "MIC"), "breakpoint_R", drop = TRUE]))))
expect_false(any(is.na(as.disk(rsi_translation[which(rsi_translation$method == "DISK"), "breakpoint_S", drop = TRUE]))))
expect_false(any(is.na(as.disk(rsi_translation[which(rsi_translation$method == "DISK"), "breakpoint_R", drop = TRUE]))))

# antibiotic names must always be coercible to their original AB code
expect_identical(as.ab(antibiotics$name), antibiotics$ab)

if (AMR:::pkg_is_available("tibble", also_load = FALSE)) {
  # there should be no diacritics (i.e. non ASCII) characters in the datasets (CRAN policy)
  datasets <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item", drop = TRUE]
  for (i in seq_len(length(datasets))) {
    dataset <- get(datasets[i], envir = asNamespace("AMR"))
    expect_identical(AMR:::dataset_UTF8_to_ASCII(dataset), dataset, info = datasets[i])
  }
}

df <- AMR:::MO_lookup
expect_true(nrow(df[which(df$prevalence == 1), , drop = FALSE]) < nrow(df[which(df$prevalence == 2), , drop = FALSE]))
expect_true(nrow(df[which(df$prevalence == 2), , drop = FALSE]) < nrow(df[which(df$prevalence == 3), , drop = FALSE]))
expect_true(all(c(
  "mo", "fullname",
  "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies",
  "rank", "ref", "lpsn", "gbif", "status", "source", "prevalence", "snomed",
  "kingdom_index", "fullname_lower", "g_species"
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
      uncategorised$species, " (", uncategorised$mo, ")"
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

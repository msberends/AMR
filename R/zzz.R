# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

.onLoad <- function(libname, pkgname) {
  assign(x = "AB_lookup",
         value = create_AB_lookup(),
         envir = asNamespace("AMR"))
  
  assign(x = "MO_lookup",
         value = create_MO_lookup(),
         envir = asNamespace("AMR"))
  
  assign(x = "MO.old_lookup",
         value = create_MO.old_lookup(),
         envir = asNamespace("AMR"))
  
  assign(x = "LANGUAGES_SUPPORTED",
         value = sort(c("en", unique(translations_file$lang))),
         envir = asNamespace("AMR"))

  # support for tibble headers (type_sum) and tibble columns content (pillar_shaft) without the need to depend on other packages
  # this was suggested by the developers of the vctrs package: 
  # https://github.com/r-lib/vctrs/blob/05968ce8e669f73213e3e894b5f4424af4f46316/R/register-s3.R
  s3_register("pillar::pillar_shaft", "ab")
  s3_register("tibble::type_sum", "ab")
  s3_register("pillar::pillar_shaft", "mo")
  s3_register("tibble::type_sum", "mo")
  s3_register("pillar::pillar_shaft", "rsi")
  s3_register("tibble::type_sum", "rsi")
  s3_register("pillar::pillar_shaft", "mic")
  s3_register("tibble::type_sum", "mic")
  s3_register("pillar::pillar_shaft", "disk")
  s3_register("tibble::type_sum", "disk")
  # support for frequency tables
  s3_register("cleaner::freq", "mo")
  s3_register("cleaner::freq", "rsi")
}

.onAttach <- function(...) {
  if (!interactive() || stats::runif(1) > 0.1 || isTRUE(as.logical(getOption("AMR_silentstart", FALSE)))) {
    return()
  }
  packageStartupMessage("Thank you for using the AMR package! ",
                        "If you have a minute, please anonymously fill in this short questionnaire to improve the package and its functionalities:",
                        "\nhttps://msberends.github.io/AMR/survey.html",
                        "\n[ prevent his notice with suppressPackageStartupMessages(library(AMR)) or use options(AMR_silentstart = TRUE) ]")
}

create_AB_lookup <- function() {
  AB_lookup <- AMR::antibiotics
  AB_lookup$generalised_name <- generalise_antibiotic_name(AB_lookup$name)
  AB_lookup$generalised_synonyms <- lapply(AB_lookup$synonyms, generalise_antibiotic_name)
  AB_lookup$generalised_abbreviations <- lapply(AB_lookup$abbreviations, generalise_antibiotic_name)
  AB_lookup$generalised_loinc <- lapply(AB_lookup$loinc, generalise_antibiotic_name)
  AB_lookup
}

create_MO_lookup <- function() {
  MO_lookup <- AMR::microorganisms
  
  MO_lookup$kingdom_index <- NA_real_
  MO_lookup[which(MO_lookup$kingdom == "Bacteria" | MO_lookup$mo == "UNKNOWN"), "kingdom_index"] <- 1
  MO_lookup[which(MO_lookup$kingdom == "Fungi"), "kingdom_index"] <- 2
  MO_lookup[which(MO_lookup$kingdom == "Protozoa"), "kingdom_index"] <- 3
  MO_lookup[which(MO_lookup$kingdom == "Archaea"), "kingdom_index"] <- 4
  # all the rest
  MO_lookup[which(is.na(MO_lookup$kingdom_index)), "kingdom_index"] <- 5
  
  # use this paste instead of `fullname` to work with Viridans Group Streptococci, etc.
  MO_lookup$fullname_lower <- tolower(trimws(paste(MO_lookup$genus, 
                                                   MO_lookup$species,
                                                   MO_lookup$subspecies)))
  ind <- MO_lookup$genus == "" | grepl("^[(]unknown ", MO_lookup$fullname)
  MO_lookup[ind, "fullname_lower"] <- tolower(MO_lookup[ind, "fullname"])
  MO_lookup$fullname_lower <- trimws(gsub("[^.a-z0-9/ \\-]+", "", MO_lookup$fullname_lower))
  
  # add a column with only "e coli" like combinations
  MO_lookup$g_species <- gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO_lookup$fullname_lower)
  
  # so arrange data on prevalence first, then kingdom, then full name
  MO_lookup[order(MO_lookup$prevalence, MO_lookup$kingdom_index, MO_lookup$fullname_lower), ]
}

create_MO.old_lookup <- function() {
  MO.old_lookup <- AMR::microorganisms.old
  MO.old_lookup$fullname_lower <- trimws(gsub("[^.a-z0-9/ \\-]+", "", tolower(trimws(MO.old_lookup$fullname))))
  
  # add a column with only "e coli"-like combinations
  MO.old_lookup$g_species <- trimws(gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO.old_lookup$fullname_lower))
  
  # so arrange data on prevalence first, then full name
  MO.old_lookup[order(MO.old_lookup$prevalence, MO.old_lookup$fullname_lower), ]
}

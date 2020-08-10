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
  assign(x = "MO_lookup",
         value = create_MO_lookup(),
         envir = asNamespace("AMR"))
  
  assign(x = "MO.old_lookup",
         value = create_MO.old_lookup(),
         envir = asNamespace("AMR"))
}

.onAttach <- function(...) {
  if (!interactive() || stats::runif(1) > 0.25 || isTRUE(as.logical(Sys.getenv("AMR_silentstart", FALSE)))) {
    return()
  }
  packageStartupMessage("Thank you for using the AMR package! ",
                        "If you have a minute, please anonymously fill in this short questionnaire to improve the package and its functionalities:",
                        "\nhttps://msberends.github.io/AMR/survey.html",
                        "\n[ permanently turn this message off with: Sys.setenv(AMR_silentstart = TRUE) ]")
}

create_MO_lookup <- function() {
  MO_lookup <- AMR::microorganisms
  
  MO_lookup$kingdom_index <- 99
  MO_lookup[which(MO_lookup$kingdom == "Bacteria" | MO_lookup$mo == "UNKNOWN"), "kingdom_index"] <- 1
  MO_lookup[which(MO_lookup$kingdom == "Fungi"), "kingdom_index"] <- 2
  MO_lookup[which(MO_lookup$kingdom == "Protozoa"), "kingdom_index"] <- 3
  MO_lookup[which(MO_lookup$kingdom == "Archaea"), "kingdom_index"] <- 4
  
  # use this paste instead of `fullname` to work with Viridans Group Streptococci, etc.
  MO_lookup$fullname_lower <- tolower(trimws(paste(MO_lookup$genus, 
                                                   MO_lookup$species,
                                                   MO_lookup$subspecies)))
  MO_lookup[MO_lookup$genus == "" | grepl("^[(]unknown ", MO_lookup$fullname), "fullname_lower"] <- tolower(trimws(MO_lookup[MO_lookup$genus == "" | grepl("^[(]unknown ", MO_lookup$fullname), 
                                                                                                                             "fullname"]))
  MO_lookup$fullname_lower <- gsub("[^.a-z0-9/ \\-]+", "", MO_lookup$fullname_lower)
  
  # add a column with only "e coli" like combinations
  MO_lookup$g_species <- gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO_lookup$fullname_lower)
  
  # so arrange data on prevalence first, then kingdom, then full name
  MO_lookup[order(MO_lookup$prevalence, MO_lookup$kingdom_index, MO_lookup$fullname_lower), ]
}

create_MO.old_lookup <- function() {
  MO.old_lookup <- AMR::microorganisms.old
  MO.old_lookup$fullname_lower <- gsub("[^.a-z0-9/ \\-]+", "", tolower(trimws(MO.old_lookup$fullname)))
  
  # add a column with only "e coli" like combinations
  MO.old_lookup$g_species <- gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO.old_lookup$fullname_lower)
  
  # so arrange data on prevalence first, then full name
  MO.old_lookup[order(MO.old_lookup$prevalence, MO.old_lookup$fullname_lower), ]
}

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

# Run this file to update the languages used in the packages:
# source("data-raw/_language_update.R")

if (!file.exists("DESCRIPTION") || !"Package: AMR" %in% readLines("DESCRIPTION")) {
  stop("Be sure to run this script in the root location of the AMR package folder.\n",
    "Working directory expected to contain the DESCRIPTION file of the AMR package.\n",
    "Current working directory: ", getwd(),
    call. = FALSE
  )
}

# save old global env to restore later
lang_env <- new.env(hash = FALSE)

# load current internal data into new env
load("R/sysdata.rda", envir = lang_env)

# replace language objects with updates
message("Reading translation file...")
lang_env$TRANSLATIONS <- utils::read.delim(
  file = "data-raw/translations.tsv",
  sep = "\t",
  stringsAsFactors = FALSE,
  header = TRUE,
  blank.lines.skip = TRUE,
  fill = TRUE,
  strip.white = TRUE,
  encoding = "UTF-8",
  fileEncoding = "UTF-8",
  na.strings = c(NA, "", NULL),
  allowEscapes = TRUE, # else "\\1" will be imported as "\\\\1"
  quote = ""
)
lang_env$TRANSLATIONS <- lang_env$TRANSLATIONS[, which(colnames(lang_env$TRANSLATIONS) != "en"), drop = FALSE]

lang_env$LANGUAGES_SUPPORTED_NAMES <- c(
  list(en = list(exonym = "English", endonym = "English")),
  lapply(
    lang_env$TRANSLATIONS[, which(nchar(colnames(lang_env$TRANSLATIONS)) == 2), drop = FALSE],
    function(x) list(exonym = x[1], endonym = x[2])
  )
)

lang_env$LANGUAGES_SUPPORTED <- names(lang_env$LANGUAGES_SUPPORTED_NAMES)

# save env to internal package data
# usethis::use_data() does not allow to save a list :(
message("Saving to internal data...")
save(
  list = names(lang_env),
  file = "R/sysdata.rda",
  ascii = FALSE,
  version = 2,
  compress = "xz",
  envir = lang_env
)

rm(lang_env)

if ("roxygen2" %in% utils::installed.packages()) {
  message("Updating package documentation...")
  suppressMessages(roxygen2::roxygenise(package.dir = "."))
} else {
  message("NOTE: please install the roxygen2 package to update package documentation, and run this script again.")
}

message("Done!")

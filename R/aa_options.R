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

#' Options for the AMR package
#'
#' This is an overview of all the package-specific [options()] you can set in the `AMR` package. 
#' @section Options:
#' * `AMR_custom_ab` \cr Allows to use custom antimicrobial drugs with this package. This is explained in [add_custom_antimicrobials()].
#' * `AMR_custom_mo` \cr Allows to use custom microorganisms with this package. This is explained in [add_custom_microorganisms()].
#' * `AMR_eucastrules` \cr Used for setting the default types of rules for [eucast_rules()] function, must be one or more of: `"breakpoints"`, `"expert"`, `"other"`, `"custom"`, `"all"`, and defaults to `c("breakpoints", "expert")`.
#' * `AMR_guideline` \cr Used for setting the default guideline for interpreting MIC values and disk diffusion diameters with [as.sir()]. Can be only the guideline name (e.g., `"CLSI"`) or the name with a year (e.g. `"CLSI 2019"`). The default is \code{"`r clinical_breakpoints$guideline[1]`"}. Supported guideline are currently EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`).
#' * `AMR_ignore_pattern` \cr A [regular expression][base::regex] to define input that must be ignored in [as.mo()] and all [`mo_*`][mo_property()] functions.
#' * `AMR_include_PKPD` \cr A [logical] to use in [as.sir()], to indicate that PK/PD clinical breakpoints must be applied as a last resort, defaults to `TRUE`.
#' * `AMR_include_screening` \cr A [logical] to use in [as.sir()], to indicate that clinical breakpoints for screening are allowed, defaults to `FALSE`.
#' * `AMR_keep_synonyms` \cr A [logical] to use in [as.mo()] and all [`mo_*`][mo_property()] functions, to indicate if old, previously valid taxonomic names must be preserved and not be corrected to currently accepted names.
#' * `AMR_locale` \cr A language to use for the `AMR` package, can be one of these supported language names or ISO-639-1 codes: `r vector_or(paste0(sapply(LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]), " (" , LANGUAGES_SUPPORTED, ")"), quotes = FALSE, sort = FALSE)`.
#' * `AMR_mo_source` \cr A file location for a manual code list to be used in [as.mo()] and all [`mo_*`][mo_property()] functions. This is explained in [set_mo_source()].
#' 
#' @section Saving Settings Between Sessions:
#' Settings in \R are not saved globally and are thus lost when \R is exited. You can save your options to your own `.Rprofile` file, which is a user-specific file. You can edit it using:
#' 
#' ```r
#'   utils::file.edit("~/.Rprofile")
#' ```
#' 
#' In this file, you can set options such as:
#' 
#' ```r
#'  options(AMR_locale = "pt")
#'  options(AMR_include_PKPD = TRUE)
#'  ```
#'  
#' to add Portuguese language support of antibiotics, and allow PK/PD rules when interpreting MIC values with [as.sir()].
#' 
#' ### Share Options Within Team
#' 
#' For a more global approach, e.g. within a data team, save an options file to a remote file location, such as a shared network drive. This would work in this way:
#' 
#' 1. Save a plain text file to e.g. "X:/team_folder/R_options.R" and fill it with preferred settings.
#' 
#' 2. For each user, open the `.Rprofile` file using `utils::file.edit("~/.Rprofile")` and put in there:
#' 
#'    ```r
#'      source("X:/team_folder/R_options.R")
#'    ```
#'    
#' 3. Reload R/RStudio and check the settings with [getOption()], e.g. `getOption("AMR_locale")` if you have set that value.
#' 
#' Now the team settings are configured in only one place, and can be maintained there.
#' @keywords internal
#' @name AMR-options
NULL

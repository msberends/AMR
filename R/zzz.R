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

# set up package environment, used by numerous AMR functions
AMR_env <- new.env(hash = TRUE, parent = emptyenv())
AMR_env$mo_uncertainties <- data.frame(
  original_input = character(0),
  input = character(0),
  fullname = character(0),
  mo = character(0),
  candidates = character(0),
  minimum_matching_score = integer(0),
  keep_synonyms = logical(0),
  stringsAsFactors = FALSE
)
AMR_env$mo_renamed <- list()
AMR_env$mo_previously_coerced <- data.frame(
  x = character(0),
  mo = character(0),
  stringsAsFactors = FALSE
)
AMR_env$ab_previously_coerced <- data.frame(
  x = character(0),
  ab = character(0),
  stringsAsFactors = FALSE
)
AMR_env$av_previously_coerced <- data.frame(
  x = character(0),
  av = character(0),
  stringsAsFactors = FALSE
)
AMR_env$sir_interpretation_history <- data.frame(
  datetime = Sys.time()[0],
  index = integer(0),
  method = character(0),
  ab_given = character(0),
  mo_given = character(0),
  host_given = character(0),
  input_given = character(0),
  ab = set_clean_class(character(0), c("ab", "character")),
  mo = set_clean_class(character(0), c("mo", "character")),
  host = character(0),
  input = character(0),
  outcome = NA_sir_[0],
  notes = character(0),
  guideline = character(0),
  ref_table = character(0),
  uti = logical(0),
  breakpoint_S_R = character(0),
  site = character(0),
  stringsAsFactors = FALSE
)


AMR_env$custom_ab_codes <- character(0)
AMR_env$custom_mo_codes <- character(0)
AMR_env$is_dark_theme <- NULL
AMR_env$supports_colour <- NULL
AMR_env$chmatch <- import_fn("chmatch", "data.table", error_on_fail = FALSE)
AMR_env$chin <- import_fn("%chin%", "data.table", error_on_fail = FALSE)

# take cli symbols and error function if available
AMR_env$bullet_icon <- import_fn("symbol", "cli", error_on_fail = FALSE)$bullet %or% "*"
AMR_env$ellipsis_icon <- import_fn("symbol", "cli", error_on_fail = FALSE)$ellipsis %or% "..."
AMR_env$info_icon <- import_fn("symbol", "cli", error_on_fail = FALSE)$info %or% "i"
AMR_env$sup_1_icon <- import_fn("symbol", "cli", error_on_fail = FALSE)$sup_1 %or% "*"
AMR_env$cli_abort <- import_fn("cli_abort", "cli", error_on_fail = FALSE)
AMR_env$cross_icon <- if (isTRUE(base::l10n_info()$`UTF-8`)) "\u00d7" else "x"

.onLoad <- function(libname, pkgname) {
  # if mo source exists, fire it up (see mo_source())
  if (tryCatch(file.exists(getOption("AMR_mo_source", "~/mo_source.rds")), error = function(e) FALSE)) {
    try(invisible(get_mo_source()), silent = TRUE)
  }
  # be sure to print tibbles as tibbles
  if (pkg_is_available("tibble")) {
    try(loadNamespace("tibble"), silent = TRUE)
  }

  # reference data - they have additional data to improve algorithm speed
  # they cannot be part of R/sysdata.rda since CRAN thinks it would make the package too large (+3 MB)
  if (NROW(AB_LOOKUP) != NROW(AMR::antimicrobials)) {
    # antimicrobials data set was updated - run create_AB_AV_lookup() again
    AB_LOOKUP <- create_AB_AV_lookup(AMR::antimicrobials)
  }

  AMR_env$AB_lookup <- cbind(AMR::antimicrobials, AB_LOOKUP)
  AMR_env$AV_lookup <- cbind(AMR::antivirals, AV_LOOKUP)
}

.onAttach <- function(libname, pkgname) {
  # if custom ab option is available, load it
  if (!is.null(getOption("AMR_custom_ab")) && file.exists(getOption("AMR_custom_ab", default = ""))) {
    if (getOption("AMR_custom_ab") %unlike% "[.]rds$") {
      packageStartupMessage("The file with custom antimicrobials must be an RDS file. Set the option `AMR_custom_ab` to another path.")
    } else {
      packageStartupMessage("Adding custom antimicrobials from '", getOption("AMR_custom_ab"), "'...", appendLF = FALSE)
      x <- readRDS_AMR(getOption("AMR_custom_ab"))
      tryCatch(
        {
          suppressWarnings(suppressMessages(add_custom_antimicrobials(x)))
          packageStartupMessage("OK.")
        },
        error = function(e) packageStartupMessage("Failed: ", conditionMessage(e))
      )
    }
  }
  # if custom mo option is available, load it
  if (!is.null(getOption("AMR_custom_mo")) && file.exists(getOption("AMR_custom_mo", default = ""))) {
    if (getOption("AMR_custom_mo") %unlike% "[.]rds$") {
      packageStartupMessage("The file with custom microorganisms must be an RDS file. Set the option `AMR_custom_mo` to another path.")
    } else {
      packageStartupMessage("Adding custom microorganisms from '", getOption("AMR_custom_mo"), "'...", appendLF = FALSE)
      x <- readRDS_AMR(getOption("AMR_custom_mo"))
      tryCatch(
        {
          suppressWarnings(suppressMessages(add_custom_microorganisms(x)))
          packageStartupMessage("OK.")
        },
        error = function(e) packageStartupMessage("Failed: ", conditionMessage(e))
      )
    }
  }
}

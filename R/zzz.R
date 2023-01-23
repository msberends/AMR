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

# set up package environment, used by numerous AMR functions
AMR_env <- new.env(hash = FALSE)
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
  ab_input = character(0),
  ab_considered = character(0),
  mo_input = character(0),
  mo_considered = character(0),
  guideline = character(0),
  ref_table = character(0),
  method = character(0),
  breakpoint_S = double(0),
  breakpoint_R = double(0),
  input = double(0),
  interpretation = character(0),
  stringsAsFactors = FALSE
)
AMR_env$custom_ab_codes <- character(0)
AMR_env$custom_mo_codes <- character(0)
AMR_env$is_dark_theme <- NULL

# determine info icon for messages
utf8_supported <- isTRUE(base::l10n_info()$`UTF-8`)
is_latex <- tryCatch(import_fn("is_latex_output", "knitr", error_on_fail = FALSE)(),
  error = function(e) FALSE
)
if (utf8_supported && !is_latex) {
  # \u2139 is a symbol officially named 'information source'
  AMR_env$info_icon <- "\u2139"
  AMR_env$bullet_icon <- "\u2022"
} else {
  AMR_env$info_icon <- "i"
  AMR_env$bullet_icon <- "*"
}

.onLoad <- function(lib, pkg) {
  # Support for tibble headers (type_sum) and tibble columns content (pillar_shaft)
  # without the need to depend on other packages. This was suggested by the
  # developers of the vctrs package:
  # https://github.com/r-lib/vctrs/blob/05968ce8e669f73213e3e894b5f4424af4f46316/R/register-s3.R
  s3_register("pillar::pillar_shaft", "ab")
  s3_register("pillar::pillar_shaft", "av")
  s3_register("pillar::pillar_shaft", "mo")
  s3_register("pillar::pillar_shaft", "sir")
  s3_register("pillar::pillar_shaft", "rsi") # remove in a later version
  s3_register("pillar::pillar_shaft", "mic")
  s3_register("pillar::pillar_shaft", "disk")
  s3_register("pillar::type_sum", "ab")
  s3_register("pillar::type_sum", "av")
  s3_register("pillar::type_sum", "mo")
  s3_register("pillar::type_sum", "sir")
  s3_register("pillar::type_sum", "rsi") # remove in a later version
  s3_register("pillar::type_sum", "mic")
  s3_register("pillar::type_sum", "disk")
  # Support for frequency tables from the cleaner package
  s3_register("cleaner::freq", "mo")
  s3_register("cleaner::freq", "sir")
  # Support for skim() from the skimr package
  if (pkg_is_available("skimr", also_load = FALSE, min_version = "2.0.0")) {
    s3_register("skimr::get_skimmers", "mo")
    s3_register("skimr::get_skimmers", "sir")
    s3_register("skimr::get_skimmers", "mic")
    s3_register("skimr::get_skimmers", "disk")
  }
  # Support for autoplot() from the ggplot2 package
  s3_register("ggplot2::autoplot", "sir")
  s3_register("ggplot2::autoplot", "mic")
  s3_register("ggplot2::autoplot", "disk")
  s3_register("ggplot2::autoplot", "resistance_predict")
  # Support for fortify from the ggplot2 package
  s3_register("ggplot2::fortify", "sir")
  s3_register("ggplot2::fortify", "mic")
  s3_register("ggplot2::fortify", "disk")
  # Support vctrs package for use in e.g. dplyr verbs
  # S3: ab_selector
  s3_register("vctrs::vec_ptype2", "character.ab_selector")
  s3_register("vctrs::vec_ptype2", "ab_selector.character")
  s3_register("vctrs::vec_cast", "character.ab_selector")
  # S3: ab_selector_any_all
  s3_register("vctrs::vec_ptype2", "logical.ab_selector_any_all")
  s3_register("vctrs::vec_ptype2", "ab_selector_any_all.logical")
  s3_register("vctrs::vec_cast", "logical.ab_selector_any_all")
  # S3: ab
  s3_register("vctrs::vec_ptype2", "character.ab")
  s3_register("vctrs::vec_ptype2", "ab.character")
  s3_register("vctrs::vec_cast", "character.ab")
  s3_register("vctrs::vec_cast", "ab.character")
  # S3: av
  s3_register("vctrs::vec_ptype2", "character.av")
  s3_register("vctrs::vec_ptype2", "av.character")
  s3_register("vctrs::vec_cast", "character.av")
  s3_register("vctrs::vec_cast", "av.character")
  # S3: mo
  s3_register("vctrs::vec_ptype2", "character.mo")
  s3_register("vctrs::vec_ptype2", "mo.character")
  s3_register("vctrs::vec_cast", "character.mo")
  s3_register("vctrs::vec_cast", "mo.character")
  # S3: disk
  s3_register("vctrs::vec_ptype2", "integer.disk")
  s3_register("vctrs::vec_ptype2", "disk.integer")
  s3_register("vctrs::vec_cast", "integer.disk")
  s3_register("vctrs::vec_cast", "disk.integer")
  s3_register("vctrs::vec_cast", "double.disk")
  s3_register("vctrs::vec_cast", "disk.double")
  s3_register("vctrs::vec_cast", "character.disk")
  s3_register("vctrs::vec_cast", "disk.character")
  # S3: mic
  s3_register("vctrs::vec_cast", "character.mic")
  s3_register("vctrs::vec_cast", "double.mic")
  s3_register("vctrs::vec_cast", "mic.character")
  s3_register("vctrs::vec_cast", "mic.double")
  s3_register("vctrs::vec_math", "mic")
  # S3: sir
  s3_register("vctrs::vec_ptype2", "character.sir")
  s3_register("vctrs::vec_ptype2", "sir.character")
  s3_register("vctrs::vec_cast", "character.sir")
  s3_register("vctrs::vec_cast", "sir.character")

  # if mo source exists, fire it up (see mo_source())
  if (tryCatch(file.exists(getOption("AMR_mo_source", "~/mo_source.rds")), error = function(e) FALSE)) {
    try(invisible(get_mo_source()), silent = TRUE)
  }
  # be sure to print tibbles as tibbles
  if (pkg_is_available("tibble", also_load = FALSE)) {
    try(loadNamespace("tibble"), silent = TRUE)
  }

  # reference data - they have additional to improve algorithm speed
  # they cannot be part of R/sysdata.rda since CRAN thinks it would make the package too large (+3 MB)
  AMR_env$AB_lookup <- cbind(AMR::antibiotics, AB_LOOKUP)
  AMR_env$AV_lookup <- cbind(AMR::antivirals, AV_LOOKUP)
}

.onAttach <- function(lib, pkg) {
  # if custom ab option is available, load it
  if (!is.null(getOption("AMR_custom_ab")) && file.exists(getOption("AMR_custom_ab", default = ""))) {
    packageStartupMessage("Adding custom antimicrobials from '", getOption("AMR_custom_ab"), "'...", appendLF = FALSE)
    x <- readRDS2(getOption("AMR_custom_ab"))
    tryCatch(
      {
        suppressWarnings(suppressMessages(add_custom_antimicrobials(x)))
        packageStartupMessage("OK.")
      },
      error = function(e) packageStartupMessage("Failed: ", e$message)
    )
  }
  # if custom mo option is available, load it
  if (!is.null(getOption("AMR_custom_mo")) && file.exists(getOption("AMR_custom_mo", default = ""))) {
    packageStartupMessage("Adding custom microorganisms from '", getOption("AMR_custom_mo"), "'...", appendLF = FALSE)
    x <- readRDS2(getOption("AMR_custom_mo"))
    tryCatch(
      {
        suppressWarnings(suppressMessages(add_custom_microorganisms(x)))
        packageStartupMessage("OK.")
      },
      error = function(e) packageStartupMessage("Failed: ", e$message)
    )
  }
}

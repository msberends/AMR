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

#' Add Manual Antimicrobials to This Package
#' 
#' With [add_custom_antimicrobials()] you can add your own manual antimicrobial codes to the `AMR` package.
#' @param x a [data.frame] resembling the [antibiotics] data set, at least containing columns "ab" and "name"
#' @details Due to how \R works, the [add_custom_antimicrobials()] function has to be run in every \R session - added antimicrobials are not stored between sessions and are thus lost when \R is exited. It is possible to save the antimicrobial additions to your `.Rprofile` file to circumvent this, for example:
#' 
#' ```r
#' library(AMR)
#' add_custom_antimicrobials(
#'   data.frame(ab = "TEST",
#'              name = "Test Antibiotic",
#'              group = "Test Group")
#' )
#' ```
#' 
#' Use [clear_custom_antimicrobials()] to clear the previously added antimicrobials.
#' @rdname add_custom_antimicrobials
#' @export
#' @examples 
#' \donttest{
#' 
#' # returns NA and throws a warning:
#' as.ab("test")
#' 
#' # now a manual entry - it will be considered by as.ab() and
#' # all ab_*() functions
#' add_custom_antimicrobials(
#'   data.frame(ab = "TEST",
#'              name = "Test Antibiotic",
#'              group = "Test Group")
#' )
#' 
#' "test" is now a new antibiotic:
#' as.ab("test")
#' ab_name("test")
#' ab_group("test")
#' }
add_custom_antimicrobials <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  stop_ifnot(all(c("ab", "name") %in% colnames(x)),
             "`x` must contain columns \"ab\" and \"name\".")
  stop_if(any(x$ab %in% AB_lookup$ab),
          "Antimicrobial code(s) ", vector_and(x$ab[x$ab %in% AB_lookup$ab]), " already exist in the internal `antibiotics` data set.")
  
  x <- x[, colnames(AB_lookup)[colnames(AB_lookup) %in% colnames(x)], drop = FALSE]
  x$generalised_name <- generalise_antibiotic_name(x$name)
  x$generalised_all <- as.list(x$generalised_name)
  
  bind_rows <- import_fn("bind_rows", "dplyr", error_on_fail = FALSE)
  if (!is.null(bind_rows)) {
    new_df <- bind_rows(AB_lookup, x)
  } else {
    new_df <- rbind(AB_lookup, x, stringsAsFactors = FALSE)
  }
  
  assignInNamespace(x = "AB_lookup",
                    value = new_df,
                    ns = asNamespace("AMR"))
  message_("Added ", nr2char(nrow(x)), " record", ifelse(nrow(x) > 1, "s", ""), " to internal `antibiotics` data set.")
}

#' @rdname add_custom_antimicrobials
#' @export
clear_custom_antimicrobials <- function() {
  assignInNamespace(x = "AB_lookup",
                    value = create_AB_lookup(),
                    ns = asNamespace("AMR"))
  message_("Manual antimicrobials cleared.")
}

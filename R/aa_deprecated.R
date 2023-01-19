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

#' Deprecated Functions
#'
#' These functions are so-called '[Deprecated]'. **They will be removed in a future release.** Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @keywords internal
#' @name AMR-deprecated
#' @rdname AMR-deprecated
#' @export
as.rsi <- function(...) {
  deprecation_warning("as.rsi", "as.sir")
  as.sir(...)
}

#' @rdname AMR-deprecated
#' @export
is.rsi.eligible <- function(...) {
  deprecation_warning("is.rsi.eligible", "is_sir_eligible")
  is_sir_eligible(...)
}

# NAMESPACE NALOPEN

# will be exported using s3_register() in R/zzz.R
pillar_shaft.rsi <- function(x, ...) {
  out <- trimws(format(x))
  if (has_colour()) {
    # colours will anyway not work when has_colour() == FALSE,
    # but then the indentation should also not be applied
    out[is.na(x)] <- font_grey("  NA")
    out[x == "S"] <- font_green_bg("  S  ")
    out[x == "I"] <- font_orange_bg("  I  ")
    out[x == "R"] <- font_red_bg("  R  ")
  }
  create_pillar_column(out, align = "left", width = 5)
}
type_sum.rsi <- function(x, ...) {
  deprecation_warning("as.rsi", "as.sir", "Transform your old 'rsi' class to the new 'sir' with `as.sir()` using e.g.:\n  your_data %>% mutate_if(~inherits(.x, \"rsi\"), as.sir)")
  "rsi"
}

#' @method print rsi
#' @export
#' @noRd
print.rsi <- function(x, ...) {
  deprecation_warning("as.rsi", "as.sir", "Transform your old 'rsi' class to the new 'sir' with `as.sir()`")
  print(x, ...)
}

deprecation_warning <- function(old, new = NULL, extra_msg = NULL) {
  env <- paste0("deprecated_", old)
  if (!env %in% names(AMR_env)) {
    AMR_env[[paste0("deprecated_", old)]] <- 1
    warning_(ifelse(is.null(new), 
                    paste0("The `", old, "()` function is no longer in use"),
                    paste0("The `", old, "()` function has been replaced with `", new, "()`")),
             ", see `?AMR-deprecated`.",
             ifelse(!is.null(extra_msg),
                    paste0(" ", extra_msg),
                    ""),
             "\nThis warning will be shown once per session.")
  }
}

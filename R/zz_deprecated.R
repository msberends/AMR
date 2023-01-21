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
NA_rsi_ <- set_clean_class(factor(NA_character_, levels = c("S", "I", "R"), ordered = TRUE),
                           new_class = c("rsi", "ordered", "factor"))
#' @rdname AMR-deprecated
#' @export
as.rsi <- function(x, ...) {
  deprecation_warning("as.rsi", "as.sir")
  UseMethod("as.rsi")
}
#' @noRd
#' @export
as.rsi.default <- function(...) {
  as.sir.default(...)
}
#' @noRd
#' @export
as.rsi.mic <- function(...) {
  as.sir.mic(...)
}
#' @noRd
#' @export
as.rsi.disk <- function(...) {
  as.sir.disk(...)
}
#' @noRd
#' @export
as.rsi.data.frame <- function(...) {
  as.sir.data.frame(...)
}

#' @rdname AMR-deprecated
#' @export
facet_rsi <- function(...) {
  deprecation_warning("facet_rsi", "facet_sir")
  facet_sir(...)
}
#' @rdname AMR-deprecated
#' @export
geom_rsi <- function(...) {
  deprecation_warning("geom_rsi", "geom_sir")
  geom_sir(...)
}
#' @rdname AMR-deprecated
#' @export
ggplot_rsi <- function(...) {
  deprecation_warning("ggplot_rsi", "ggplot_sir")
  ggplot_sir(...)
}
#' @rdname AMR-deprecated
#' @export
ggplot_rsi_predict <- function(...) {
  deprecation_warning("ggplot_rsi_predict", "ggplot_sir_predict")
  ggplot_sir_predict(...)
}
#' @rdname AMR-deprecated
#' @export
is.rsi <- function(x, ...) {
  # this is an exception, so mutate_if(is.rsi, as.sir) can be used
  if (inherits(x, "data.frame")) {
    unname(vapply(FUN.VALUE = logical(1), x, is.rsi))
  } else {
    inherits(x, "rsi")
  }
}
#' @rdname AMR-deprecated
#' @export
is.rsi.eligible <- function(...) {
  deprecation_warning("is.rsi.eligible", "is_sir_eligible")
  is_sir_eligible(...)
}
#' @rdname AMR-deprecated
#' @export
labels_rsi_count <- function(...) {
  deprecation_warning("labels_rsi_count", "labels_sir_count")
  labels_sir_count(...)
}
#' @rdname AMR-deprecated
#' @export
n_rsi <- function(...) {
  deprecation_warning("n_rsi", "n_sir")
  n_sir(...)
}
#' @rdname AMR-deprecated
#' @export
random_rsi <- function(...) {
  deprecation_warning("random_rsi", "random_sir")
  random_sir(...)
}
#' @rdname AMR-deprecated
#' @export
rsi_df <- function(...) {
  deprecation_warning("rsi_df", "sir_df")
  sir_df(...)
}
#' @rdname AMR-deprecated
#' @export
rsi_predict <- function(...) {
  deprecation_warning("rsi_predict", "sir_predict")
  sir_predict(...)
}
#' @rdname AMR-deprecated
#' @export
scale_rsi_colours <- function(...) {
  deprecation_warning("scale_rsi_colours", "scale_sir_colours")
  scale_sir_colours(...)
}
#' @rdname AMR-deprecated
#' @export
theme_rsi <- function(...) {
  deprecation_warning("theme_rsi", "theme_sir")
  theme_sir(...)
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.rsi <- pillar_shaft.sir
type_sum.rsi <- function(x, ...) {
  deprecation_warning(extra_msg = "* Transform your old 'rsi' class to the new 'sir' class with `as.sir()` using e.g.:\n  your_data %>% mutate_if(is.rsi, as.sir)")
  paste0("rsi", font_bold(font_red("[!]")))
}

#' @method print rsi
#' @export
#' @noRd
print.rsi <- function(x, ...) {
  deprecation_warning(extra_msg = "Transform your old 'rsi' class to the new 'sir' class with `as.sir()`")
  cat("Class 'rsi'", font_bold(font_red("[!]\n")))
  print(as.character(x), quote = FALSE)
}

#' @noRd
#' @export
`[<-.rsi` <- `[<-.sir`
#' @noRd
#' @export
`[[<-.rsi` <- `[[<-.sir`
#' @noRd
#' @export
barplot.rsi <- barplot.sir
#' @noRd
#' @export
c.rsi <- c.sir
#' @noRd
#' @export
droplevels.rsi <- droplevels.sir
#' @noRd
#' @export
plot.rsi <- plot.sir
#' @noRd
#' @export
rep.rsi <- rep.sir
#' @noRd
#' @export
summary.rsi <- summary.sir
#' @noRd
#' @export
unique.rsi <- unique.sir

deprecation_warning <- function(old = NULL, new = NULL, extra_msg = NULL) {
  if (is.null(old)) {
    warning_(extra_msg)
  } else {
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
}

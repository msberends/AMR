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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

# These are all S3 implementations for the vctrs package,
# that is used internally by tidyverse packages such as dplyr.
# They are to convert AMR-specific classes to bare characters and integers.
# All of them will be exported using s3_register() in R/zzz.R when loading the package.

# see https://github.com/tidyverse/dplyr/issues/5955 why this is required

# S3: amr_selector ----
# this does not need a .default method since it's used internally only
vec_ptype2.character.amr_selector <- function(x, y, ...) {
  x
}
vec_ptype2.amr_selector.character <- function(x, y, ...) {
  y
}
vec_cast.character.amr_selector <- function(x, to, ...) {
  unclass(x)
}

# S3: amr_selector_any_all ----
# this does not need a .default method since it's used internally only
vec_ptype2.logical.amr_selector_any_all <- function(x, y, ...) {
  x
}
vec_ptype2.amr_selector_any_all.logical <- function(x, y, ...) {
  y
}
vec_cast.logical.amr_selector_any_all <- function(x, to, ...) {
  unclass(x)
}

# S3: ab ----
vec_ptype2.ab.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  x
}
vec_ptype2.ab.ab <- function(x, y, ...) {
  x
}
vec_cast.character.ab <- function(x, to, ...) {
  as.character(x)
}
vec_cast.ab.character <- function(x, to, ...) {
  return_after_integrity_check(x, "antimicrobial drug code", as.character(AMR_env$AB_lookup$ab))
}

# S3: av ----
vec_ptype2.av.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  x
}
vec_ptype2.av.av <- function(x, y, ...) {
  x
}
vec_cast.character.av <- function(x, to, ...) {
  as.character(x)
}
vec_cast.av.character <- function(x, to, ...) {
  return_after_integrity_check(x, "antiviral drug code", as.character(AMR_env$AV_lookup$av))
}

# S3: mo ----
vec_ptype2.mo.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  x
}
vec_ptype2.mo.mo <- function(x, y, ...) {
  x
}
vec_cast.character.mo <- function(x, to, ...) {
  as.character(x)
}
vec_cast.mo.character <- function(x, to, ...) {
  add_MO_lookup_to_AMR_env()
  return_after_integrity_check(x, "microorganism code", as.character(AMR_env$MO_lookup$mo))
}

# S3: disk ----
vec_ptype_full.disk <- function(x, ...) {
  "disk"
}
vec_ptype_abbr.disk <- function(x, ...) {
  "dsk"
}
vec_ptype2.disk.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  NA_disk_[0]
}
vec_ptype2.disk.disk <- function(x, y, ...) {
  NA_disk_[0]
}
vec_cast.disk.disk <- function(x, to, ...) {
  as.disk(x)
}
vec_cast.integer.disk <- function(x, to, ...) {
  unclass(x)
}
vec_cast.disk.integer <- function(x, to, ...) {
  as.disk(x)
}
vec_cast.double.disk <- function(x, to, ...) {
  unclass(x)
}
vec_cast.disk.double <- function(x, to, ...) {
  as.disk(x)
}
vec_cast.character.disk <- function(x, to, ...) {
  unclass(x)
}
vec_cast.disk.character <- function(x, to, ...) {
  as.disk(x)
}

# S3: mic ----
vec_ptype2.mic.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  # this will make sure that currently implemented MIC levels are returned
  NA_mic_[0]
}
vec_ptype2.mic.mic <- function(x, y, ...) {
  # this will make sure that currently implemented MIC levels are returned
  NA_mic_[0]
}
vec_cast.mic.mic <- function(x, to, ...) {
  # this will make sure that currently implemented MIC levels are returned
  as.mic(x)
}
vec_cast.character.mic <- function(x, to, ...) {
  as.character(x)
}
vec_cast.double.mic <- function(x, to, ...) {
  as.double(x)
}
vec_cast.integer.mic <- function(x, to, ...) {
  as.integer(x)
}
vec_cast.factor.mic <- function(x, to, ...) {
  factor(as.character(x))
}
vec_cast.mic.double <- function(x, to, ...) {
  as.mic(x)
}
vec_cast.mic.character <- function(x, to, ...) {
  as.mic(x)
}
vec_cast.mic.integer <- function(x, to, ...) {
  as.mic(x)
}
vec_cast.mic.factor <- function(x, to, ...) {
  as.mic(x)
}
vec_math.mic <- function(.fn, x, ...) {
  .fn(as.double(x), ...)
}
vec_arith.mic <- function(op, x, y, ...) {
  vctrs::vec_arith(op, as.double(x), as.double(y))
}

# S3: sir ----
vec_ptype2.sir.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  NA_sir_[0]
}
vec_ptype2.sir.sir <- function(x, y, ...) {
  NA_sir_[0]
}
vec_ptype2.character.sir <- function(x, y, ...) {
  NA_sir_[0]
}
vec_cast.sir.sir <- function(x, to, ...) {
  # this makes sure that old SIR objects (with S/I/R) are converted to the current structure (S/SDD/I/R/NI)
  as.sir(x)
}
vec_cast.character.sir <- function(x, to, ...) {
  as.character(x)
}
vec_cast.sir.character <- function(x, to, ...) {
  as.sir(x)
}

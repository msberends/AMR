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

# These are all S3 implementations for the vctrs package,
# that is used internally by tidyverse packages such as dplyr.
# They are to convert AMR-specific classes to bare characters and integers.
# All of them will be exported using s3_register() in R/zzz.R when loading the package.

# see https://github.com/tidyverse/dplyr/issues/5955 why this is required

# S3: ab_selector
vec_ptype2.character.ab_selector <- function(x, y, ...) {
  x
}
vec_ptype2.ab_selector.character <- function(x, y, ...) {
  y
}
vec_cast.character.ab_selector <- function(x, to, ...) {
  unclass(x)
}

# S3: ab_selector_any_all
vec_ptype2.logical.ab_selector_any_all <- function(x, y, ...) {
  x
}
vec_ptype2.ab_selector_any_all.logical <- function(x, y, ...) {
  y
}
vec_cast.logical.ab_selector_any_all <- function(x, to, ...) {
  unclass(x)
}

# S3: ab
vec_ptype2.character.ab <- function(x, y, ...) {
  x
}
vec_ptype2.ab.character <- function(x, y, ...) {
  y
}
vec_cast.character.ab <- function(x, to, ...) {
  unclass(x)
}

# S3: mo
vec_ptype2.character.mo <- function(x, y, ...) {
  x
}
vec_ptype2.mo.character <- function(x, y, ...) {
  y
}
vec_cast.character.mo <- function(x, to, ...) {
  unclass(x)
}

# S3: disk
vec_ptype2.integer.disk <- function(x, y, ...) {
  x
}
vec_ptype2.disk.integer <- function(x, y, ...) {
  y
}
vec_cast.integer.disk <- function(x, to, ...) {
  unclass(x)
}

# S3: mic
vec_cast.character.mic <- function(x, to, ...) {
  as.character(x)
}
vec_cast.double.mic <- function(x, to, ...) {
  # this calls as.double.mic()
  as.double(x)
}
vec_math.mic <- function(.fn, x, ...) {
  .fn(as.double(x), ...)
}

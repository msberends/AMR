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

#' @rdname proportion
#' @export
rsi_df <- function(data,
                   translate_ab = "name",
                   language = get_AMR_locale(),
                   minimum = 30,
                   as_percent = FALSE,
                   combine_SI = TRUE,
                   confidence_level = 0.95) {
  tryCatch(
    rsi_calc_df(
      type = "both",
      data = data,
      translate_ab = translate_ab,
      language = language,
      minimum = minimum,
      as_percent = as_percent,
      combine_SI = combine_SI,
      confidence_level = confidence_level
    ),
    error = function(e) stop_(gsub("in rsi_calc_df(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

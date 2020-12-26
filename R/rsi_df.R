# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
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
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
# ==================================================================== #

#' @rdname proportion 
#' @export
rsi_df <- function(data,
                   translate_ab = "name",
                   language = get_locale(),
                   minimum = 30,
                   as_percent = FALSE,
                   combine_SI = TRUE,
                   combine_IR = FALSE) {
  rsi_calc_df(type = "both",
              data = data,
              translate_ab = translate_ab,
              language = language,
              minimum = minimum,
              as_percent = as_percent,
              combine_SI = combine_SI,
              combine_IR = combine_IR,
              combine_SI_missing = missing(combine_SI))
  
}

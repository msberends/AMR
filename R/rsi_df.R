# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' @rdname proportion
#' @rdname count
#' @export
rsi_df <- function(data,
                   translate_ab = "name",
                   language = get_locale(),
                   minimum = 30,
                   as_percent = FALSE,
                   combine_SI = TRUE,
                   combine_IR = FALSE) {
  
  proportions <- rsi_calc_df(type = "proportion",
                             data = data,
                             translate_ab = translate_ab,
                             language = language,
                             minimum = minimum,
                             as_percent = as_percent,
                             combine_SI = combine_SI,
                             combine_IR = combine_IR,
                             combine_SI_missing = missing(combine_SI))
  
  counts <- rsi_calc_df(type = "count",
                        data = data,
                        translate_ab = FALSE,
                        language = "en",
                        minimum = minimum,
                        as_percent = as_percent,
                        combine_SI = combine_SI,
                        combine_IR = combine_IR,
                        combine_SI_missing = missing(combine_SI))
  
  data.frame(proportions,
             isolates = counts$value,
             stringsAsFactors = FALSE)
  
}

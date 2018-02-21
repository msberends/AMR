# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

# No export, no Rd
"%like%" <- function(vector, pattern) {
  # Source: https://github.com/Rdatatable/data.table/blob/master/R/like.R
  if (is.factor(vector)) {
    as.integer(vector) %in% grep(pattern, levels(vector))
  } else {
    grepl(pattern, vector)
  }
}

percent <- function(x, round = 1, ...) {
  base::paste0(base::round(x * 100, digits = round), "%")
}

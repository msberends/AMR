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
  # But made it case insensitive
  if (is.factor(vector)) {
    as.integer(vector) %in% grep(pattern, levels(vector), ignore.case = TRUE)
  } else {
    grepl(pattern, vector, ignore.case = TRUE)
  }
}

# No export, no Rd
percent <- function(x, round = 1, ...) {
  base::paste0(base::round(x * 100, digits = round), "%")
}

# No export, no Rd
quasiquotate <- function(deparsed, parsed) {
  # when text: remove first and last "
  if (any(deparsed %like% '^".+"$' | deparsed %like% "^'.+'$")) {
    deparsed <- deparsed %>% substr(2, nchar(.) - 1)
  }
  # apply if needed
  if (any(!deparsed %like% '[[$:()]'
      & !deparsed %in% c('""', "''", "", # empty text
                         ".", ".data", # dplyr references
                         "TRUE", "FALSE", # logicals
                         "NA", "NaN", "NULL", # empty values
                         ls(.GlobalEnv)))) {
    deparsed
  } else {
    parsed
  }
}

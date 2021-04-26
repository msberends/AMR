# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Deprecated Functions
#'
#' These functions are so-called '[Deprecated]'. They will be removed in a future release. Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @inheritSection lifecycle Retired Lifecycle
#' @inheritSection AMR Read more on Our Website!
#' @keywords internal
#' @name AMR-deprecated
#' @export
p_symbol <- function(p, emptychar = " ") {
  .Deprecated(package = "AMR", new = "cleaner::p_symbol")
  
  p <- as.double(p)
  s <- rep(NA_character_, length(p))
  
  s[p <= 1] <- emptychar
  s[p <= 0.100] <- "."
  s[p <= 0.050] <- "*"
  s[p <= 0.010] <- "**"
  s[p <= 0.001] <- "***"
  
  s
}

#' @name AMR-deprecated
#' @export
key_antibiotics <- function(x = NULL,
                            col_mo = NULL,
                            universal_1 = guess_ab_col(x, "amoxicillin"),
                            universal_2 = guess_ab_col(x, "amoxicillin/clavulanic acid"),
                            universal_3 = guess_ab_col(x, "cefuroxime"),
                            universal_4 = guess_ab_col(x, "piperacillin/tazobactam"),
                            universal_5 = guess_ab_col(x, "ciprofloxacin"),
                            universal_6 = guess_ab_col(x, "trimethoprim/sulfamethoxazole"),
                            GramPos_1 = guess_ab_col(x, "vancomycin"),
                            GramPos_2 = guess_ab_col(x, "teicoplanin"),
                            GramPos_3 = guess_ab_col(x, "tetracycline"),
                            GramPos_4 = guess_ab_col(x, "erythromycin"),
                            GramPos_5 = guess_ab_col(x, "oxacillin"),
                            GramPos_6 = guess_ab_col(x, "rifampin"),
                            GramNeg_1 = guess_ab_col(x, "gentamicin"),
                            GramNeg_2 = guess_ab_col(x, "tobramycin"),
                            GramNeg_3 = guess_ab_col(x, "colistin"),
                            GramNeg_4 = guess_ab_col(x, "cefotaxime"),
                            GramNeg_5 = guess_ab_col(x, "ceftazidime"),
                            GramNeg_6 = guess_ab_col(x, "meropenem"),
                            warnings = TRUE,
                            ...) {
  
  .Deprecated(old = "key_antibiotics()",
              new = "key_antimicrobials()",
              package = "AMR")
  
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }

  key_antimicrobials(x = x, 
                     col_mo = col_mo, 
                     universal = c(universal_1, universal_2, universal_3, universal_4, universal_5, universal_6),
                     gram_negative = c(GramNeg_1, GramNeg_2, GramNeg_3, GramNeg_4, GramNeg_5, GramNeg_6),
                     gram_positive = c(GramPos_1, GramPos_2, GramPos_3, GramPos_4, GramPos_5, GramPos_6),
                     antifungal = NULL,
                     only_rsi_columns = FALSE,
                     ...)
}

#' @name AMR-deprecated
#' @export
key_antibiotics_equal <- function(y,
                                  z,
                                  type = "keyantimicrobials",
                                  ignore_I = TRUE,
                                  points_threshold = 2,
                                  info = FALSE,
                                  na.rm = TRUE,
                                  ...) {
  
  .Deprecated(old = "key_antibiotics_equal()",
              new = "antimicrobials_equal()",
              package = "AMR")
  
  antimicrobials_equal(y = y,
                       z = z,
                       type = type,
                       ignore_I = ignore_I,
                       points_threshold = points_threshold,
                       info = info)
}

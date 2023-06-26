# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

library(AMR)
library(dplyr)
int_resis <- data.frame(mo = microorganisms$mo, stringsAsFactors = FALSE)
for (i in seq_len(nrow(antibiotics))) {
  int_resis$new <- as.sir("S")
  colnames(int_resis)[ncol(int_resis)] <- antibiotics$ab[i]
}

int_resis <- eucast_rules(int_resis,
  eucast_rules_df = subset(
    AMR:::EUCAST_RULES_DF,
    is.na(have_these_values) & reference.version == 3.3
  ),
  info = FALSE
)

int_resis2 <- int_resis[, sapply(int_resis, function(x) any(!is.sir(x) | x == "R")), drop = FALSE] %>%
  tidyr::pivot_longer(-mo) %>%
  filter(value == "R") %>%
  select(mo, ab = name)

# remove lab drugs
untreatable <- antibiotics[which(antibiotics$name %like% "-high|EDTA|polysorbate|macromethod|screening|/nacubactam"), "ab", drop = TRUE]
# takes ages with filter()..., weird
int_resis3 <- int_resis2[which(!int_resis2$ab %in% untreatable), ]
class(int_resis3$ab) <- c("ab", "character")
int_resis3

all(int_resis3$mo %in% microorganisms$mo)
all(int_resis3$ab %in% antibiotics$ab)

intrinsic_resistant <- df_remove_nonASCII(int_resis3)
usethis::use_data(intrinsic_resistant, internal = FALSE, overwrite = TRUE, version = 2, compress = "xz")
rm(intrinsic_resistant)

# AFTER THIS:
# DO NOT FORGET TO UPDATE THE VERSION NUMBER IN mo_is_intrinsic_resistant() AND R/data.R

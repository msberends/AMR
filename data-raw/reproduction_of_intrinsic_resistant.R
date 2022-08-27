# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

library(AMR)
library(dplyr)
int_resis <- data.frame(mo = microorganisms$mo, stringsAsFactors = FALSE)
for (i in seq_len(nrow(antibiotics))) {
  int_resis$new <- as.rsi("S")
  colnames(int_resis)[ncol(int_resis)] <- antibiotics$ab[i]
}

int_resis <- eucast_rules(int_resis,
                          eucast_rules_df = subset(AMR:::EUCAST_RULES_DF,
                                                   is.na(have_these_values) & reference.version == 3.3),
                          info = FALSE)

int_resis2 <- int_resis[, sapply(int_resis, function(x) any(!is.rsi(x) | x == "R")), drop = FALSE] %>% 
  tidyr::pivot_longer(-mo) %>%
  filter(value == "R") %>% 
  select(mo, ab = name)

# remove lab drugs
untreatable <- antibiotics[which(antibiotics$name %like% "-high|EDTA|polysorbate|macromethod|screening|/nacubactam"), "ab", drop = TRUE]
int_resis2 <- int_resis2 %>% 
  filter(!ab %in% untreatable) %>% 
  arrange(mo, ab)

intrinsic_resistant <- as.data.frame(int_resis2, stringsAsFactors = FALSE)
usethis::use_data(intrinsic_resistant, internal = FALSE, overwrite = TRUE, version = 2, compress = "xz")
rm(intrinsic_resistant)

# AFTER THIS:
# DO NOT FORGET TO UPDATE THE VERSION NUMBER IN mo_is_intrinsic_resistant()

# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

library(AMR)
int_resis <- data.frame(microorganism = microorganisms$mo, stringsAsFactors = FALSE)
for (i in seq_len(nrow(antibiotics))) {
  int_resis$new <- as.rsi("S")
  colnames(int_resis)[ncol(int_resis)] <- antibiotics$name[i]
}

int_resis <- eucast_rules(int_resis, 
                          eucast_rules_df = subset(AMR:::eucast_rules_file, is.na(have_these_values)))

int_resis <- int_resis[, sapply(int_resis, function(x) any(!is.rsi(x) | x == "R"))] %>% 
  tidyr::pivot_longer(-microorganism) %>%
  filter(value == "R") %>% 
  select(microorganism, antibiotic = name)

int_resis$microorganism <- mo_name(int_resis$microorganism, language = NULL)

intrinsic_resistant <- as.data.frame(int_resis, stringsAsFactors = FALSE)
usethis::use_data(intrinsic_resistant, internal = FALSE, overwrite = TRUE, version = 2)
rm(intrinsic_resistant)

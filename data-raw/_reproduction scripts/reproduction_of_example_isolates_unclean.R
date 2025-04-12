# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

patients <- unlist(lapply(LETTERS, paste0, 1:10))

patients_table <- data.frame(
  patient_id = patients,
  gender = c(
    rep("M", 135),
    rep("F", 125)
  )
)

dates <- seq(as.Date("2011-01-01"), as.Date("2020-01-01"), by = "day")

bacteria_a <- c(
  "E. coli", "S. aureus",
  "S. pneumoniae", "K. pneumoniae"
)

bacteria_b <- c("esccol", "staaur", "strpne", "klepne")

bacteria_c <- c(
  "Escherichia coli", "Staphylococcus aureus",
  "Streptococcus pneumoniae", "Klebsiella pneumoniae"
)

ab_interpretations <- c("S", "I", "R")

ab_interpretations_messy <- c("R", "< 0.5 S", "I")

sample_size <- 1000

data_a <- data.frame(
  date = sample(dates, size = sample_size, replace = TRUE),
  hospital = "A",
  bacteria = sample(bacteria_a,
    size = sample_size, replace = TRUE,
    prob = c(0.50, 0.25, 0.15, 0.10)
  ),
  AMX = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.60, 0.05, 0.35)
  ),
  AMC = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.75, 0.10, 0.15)
  ),
  CIP = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.80, 0.00, 0.20)
  ),
  GEN = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.92, 0.00, 0.08)
  )
)


data_b <- data.frame(
  date = sample(dates, size = sample_size, replace = TRUE),
  hospital = "B",
  bacteria = sample(bacteria_b,
    size = sample_size, replace = TRUE,
    prob = c(0.50, 0.25, 0.15, 0.10)
  ),
  AMX = sample(ab_interpretations_messy,
    size = sample_size, replace = TRUE,
    prob = c(0.60, 0.05, 0.35)
  ),
  AMC = sample(ab_interpretations_messy,
    size = sample_size, replace = TRUE,
    prob = c(0.75, 0.10, 0.15)
  ),
  CIP = sample(ab_interpretations_messy,
    size = sample_size, replace = TRUE,
    prob = c(0.80, 0.00, 0.20)
  ),
  GEN = sample(ab_interpretations_messy,
    size = sample_size, replace = TRUE,
    prob = c(0.92, 0.00, 0.08)
  )
)

data_c <- data.frame(
  date = sample(dates, size = sample_size, replace = TRUE),
  hospital = "C",
  bacteria = sample(bacteria_c,
    size = sample_size, replace = TRUE,
    prob = c(0.50, 0.25, 0.15, 0.10)
  ),
  AMX = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.60, 0.05, 0.35)
  ),
  AMC = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.75, 0.10, 0.15)
  ),
  CIP = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.80, 0.00, 0.20)
  ),
  GEN = sample(ab_interpretations,
    size = sample_size, replace = TRUE,
    prob = c(0.92, 0.00, 0.08)
  )
)


example_isolates_unclean <- data_a %>%
  bind_rows(data_b, data_c)

example_isolates_unclean$patient_id <- sample(patients, size = nrow(example_isolates_unclean), replace = TRUE)

example_isolates_unclean <- example_isolates_unclean %>%
  select(patient_id, hospital, date, bacteria, everything()) %>%
  dataset_UTF8_to_ASCII()

usethis::use_data(example_isolates_unclean, overwrite = TRUE, internal = FALSE, version = 2, compress = "xz")

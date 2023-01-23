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

library(dplyr)
library(readxl)
library(cleaner)

# URL:
# https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/Dosages_v_11.0_Breakpoint_Tables.pdf
# download the PDF file, open in Adobe Acrobat and export as Excel workbook
breakpoints_version <- 12

dosage_source <- read_excel("data-raw/Dosages_v_12.0_Breakpoint_Tables.xlsx", skip = 4, na = "None") %>%
  format_names(snake_case = TRUE, penicillins = "drug") %>%
  filter(!tolower(standard_dosage) %in% c("standard dosage", "standard dosage_source", "under review")) %>%
  filter(!is.na(standard_dosage)) %>%
  # keep only one drug in the table
  arrange(desc(drug)) %>%
  mutate(drug = gsub("(.*) ([(]|iv|oral).*", "\\1", drug)) %>%
  # distinct(drug, .keep_all = TRUE) %>%
  arrange(drug) %>%
  mutate(
    ab = as.ab(drug),
    ab_name = ab_name(ab, language = NULL)
  )

dosage_source <- bind_rows(
  # oral
  dosage_source %>%
    filter(standard_dosage %like% " oral") %>%
    mutate(
      standard_dosage = gsub("oral.*", "oral", standard_dosage),
      high_dosage = if_else(high_dosage %like% "oral",
        gsub("oral.*", "oral", high_dosage),
        NA_character_
      )
    ),
  # iv
  dosage_source %>%
    filter(standard_dosage %like% " iv") %>%
    mutate(
      standard_dosage = gsub(".* or ", "", standard_dosage),
      high_dosage = if_else(high_dosage %like% "( or | iv)",
        gsub(".* or ", "", high_dosage),
        NA_character_
      )
    ),
  # im
  dosage_source %>%
    filter(standard_dosage %like% " im")
) %>%
  arrange(drug)


get_dosage_lst <- function(col_data) {
  standard <- col_data %>%
    # remove new lines
    gsub(" ?(\n|\t)+ ?", " ", .) %>%
    # keep only the first suggestion, replace all after 'or' and more informative texts
    gsub("(.*?) (or|with|loading|depending|over|by) .*", "\\1", .) %>%
    # remove (1 MU)
    gsub(" [(][0-9] [A-Z]+[)]", "", .) %>%
    # remove parentheses
    gsub("[)(]", "", .) %>%
    # remove drug names
    gsub(" [a-z]{5,99}( |$)", " ", .) %>%
    gsub(" [a-z]{5,99}( |$)", " ", .) %>%
    gsub(" (acid|dose)", "", .) # %>%
  # keep lowest value only (25-30 mg -> 25 mg)
  # gsub("[-].*? ", " ", .)

  dosage_lst <- lapply(
    strsplit(standard, " x "),
    function(x) {
      dose <- x[1]
      if (dose %like% "under") {
        dose <- NA_character_
      }
      admin <- x[2]

      list(
        dose = trimws(dose),
        dose_times = gsub("^([0-9.]+).*", "\\1", admin),
        administration = clean_character(admin),
        notes = "",
        original_txt = ""
      )
    }
  )
  for (i in seq_len(length(col_data))) {
    dosage_lst[[i]]$original_txt <- gsub("\n", " ", col_data[i])
    if (col_data[i] %like% " (or|with|loading|depending|over) ") {
      dosage_lst[[i]]$notes <- gsub("\n", " ", gsub(".* ((or|with|loading|depending|over) .*)", "\\1", col_data[i]))
    }
  }
  dosage_lst
}

standard <- get_dosage_lst(dosage_source$standard_dosage)
high <- get_dosage_lst(dosage_source$high_dosage)
uti <- get_dosage_lst(dosage_source$uncomplicated_uti)
dosage_new <- bind_rows(
  # standard dose
  data.frame(
    ab = dosage_source$ab,
    name = dosage_source$ab_name,
    type = "standard_dosage",
    dose = sapply(standard, function(x) x$dose),
    dose_times = sapply(standard, function(x) x$dose_times),
    administration = sapply(standard, function(x) x$administration),
    notes = sapply(standard, function(x) x$notes),
    original_txt = sapply(standard, function(x) x$original_txt),
    stringsAsFactors = FALSE
  ),
  # high dose
  data.frame(
    ab = dosage_source$ab,
    name = dosage_source$ab_name,
    type = "high_dosage",
    dose = sapply(high, function(x) x$dose),
    dose_times = sapply(high, function(x) x$dose_times),
    administration = sapply(high, function(x) x$administration),
    notes = sapply(high, function(x) x$notes),
    original_txt = sapply(high, function(x) x$original_txt),
    stringsAsFactors = FALSE
  ),
  # UTIs
  data.frame(
    ab = dosage_source$ab,
    name = dosage_source$ab_name,
    type = "uncomplicated_uti",
    dose = sapply(uti, function(x) x$dose),
    dose_times = sapply(uti, function(x) x$dose_times),
    administration = sapply(uti, function(x) x$administration),
    notes = sapply(uti, function(x) x$notes),
    original_txt = sapply(uti, function(x) x$original_txt),
    stringsAsFactors = FALSE
  )
) %>%
  mutate(
    eucast_version = breakpoints_version,
    dose_times = as.integer(dose_times),
    administration = gsub("([a-z]+) .*", "\\1", administration)
  ) %>%
  arrange(name, administration, type) %>%
  filter(!is.na(dose), dose != ".") %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(dosage_new) <- NULL

dosage <- bind_rows(dosage_new, AMR::dosage) %>%
  dataset_UTF8_to_ASCII()

usethis::use_data(dosage, internal = FALSE, overwrite = TRUE, version = 2)

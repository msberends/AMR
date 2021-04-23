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

library(openxlsx)
library(dplyr)
library(tidyr)
library(cleaner)
library(AMR)

# USE THIS FUNCTION TO READ THE EUCAST EXCEL FILE THAT CONTAINS THE BREAKPOINT TABLES

read_EUCAST <- function(sheet, file, guideline_name) { 
  
  message("\nGetting sheet: ", sheet)
  sheet.bak <- sheet
  
  uncertainties <- NULL
  add_uncertainties <- function(old, new) {
    if (is.null(old)) {
      new
    } else {
      bind_rows(old, new)
    }
  }
  
  raw_data <- read.xlsx(xlsxFile = file,
                        sheet = sheet,
                        colNames = FALSE,
                        skipEmptyRows = FALSE,
                        skipEmptyCols = FALSE, 
                        fillMergedCells = TRUE,
                        na.strings = c("", "-", "NA", "IE", "IP"))
  probable_rows <- suppressWarnings(raw_data %>% mutate_all(as.double) %>% summarise_all(~sum(!is.na(.))) %>% unlist() %>% max())
  if (probable_rows == 0) {
    message("NO ROWS FOUND")
    message("------------------------")
    return(NULL)
  }
  
  # in the info header in the Excel file, EUCAST mentions which genera are targeted
  if (sheet %like% "anaerob.*Gram.*posi") {
    sheet <- paste0(c("Actinomyces", "Bifidobacterium", "Clostridioides", 
                      "Clostridium", "Cutibacterium", "Eggerthella", 
                      "Eubacterium", "Lactobacillus", "Propionibacterium", 
                      "Staphylococcus saccharolyticus"), 
                    collapse = "_")
  } else if (sheet %like% "anaerob.*Gram.*nega") {
    sheet <- paste0(c("Bacteroides",
                      "Bilophila",
                      "Fusobacterium",
                      "Mobiluncus",
                      "Parabacteroides",
                      "Porphyromonas",
                      "Prevotella"), 
                    collapse = "_")
  } else if (sheet == "Streptococcus A,B,C,G") {
    sheet <- paste0(microorganisms %>%
                      filter(genus == "Streptococcus") %>%
                      mutate(lancefield = mo_name(mo, Lancefield = TRUE)) %>%
                      filter(lancefield %like% "^Streptococcus group") %>% 
                      pull(fullname), 
                    collapse = "_")
  } else if (sheet %like% "PK.*PD") {
    sheet <- "UNKNOWN"
  }
  mo_sheet <- paste0(suppressMessages(as.mo(unlist(strsplit(sheet, "_")))), collapse = "|")
  if (!is.null(mo_uncertainties())) uncertainties <- add_uncertainties(uncertainties, mo_uncertainties())
  
  set_columns_names <- function(x, cols) {
    colnames(x) <- cols[1:length(colnames(x))]
    x
  }
  
  get_mo <- function(x) {
    for (i in seq_len(length(x))) {
      y <- trimws(unlist(strsplit(x[i], "(,|and)")))
      y <- trimws(gsub("[(].*[)]", "", y))
      y <- suppressWarnings(suppressMessages(as.mo(y, allow_uncertain = FALSE)))
      if (!is.null(mo_uncertainties())) uncertainties <<- add_uncertainties(uncertainties, mo_uncertainties())
      y <- y[!is.na(y) & y != "UNKNOWN"]
      x[i] <- paste(y, collapse = "|")
    }
    x
  }
  
  MICs_with_trailing_superscript <- c(seq(from = 0.0011, to = 0.0019, by = 0.0001),
                                      seq(from = 0.031, to = 0.039, by = 0.001),
                                      seq(from = 0.061, to = 0.069, by = 0.001),
                                      seq(from = 0.1251, to = 0.1259, by = 0.0001),
                                      seq(from = 0.251, to = 0.259, by = 0.001),
                                      seq(from = 0.51, to = 0.59, by = 0.01),
                                      seq(from = 11, to = 19, by = 1),
                                      seq(from = 161, to = 169, by = 01),
                                      seq(from = 21, to = 29, by = 1),
                                      seq(from = 321, to = 329, by = 1),
                                      seq(from = 41, to = 49, by = 1),
                                      seq(from = 81, to = 89, by = 1))
  has_superscript <- function(x) {
    # because due to floating point error 0.1252 is not in: 
    # seq(from = 0.1251, to = 0.1259, by = 0.0001)
    sapply(x, function(x) any(near(x, MICs_with_trailing_superscript)))
  }
  
  has_zone_diameters <- rep(any(unlist(raw_data) %like% "zone diameter"), nrow(raw_data))

  cleaned <- raw_data %>% 
    as_tibble() %>% 
    set_columns_names(LETTERS) %>% 
    transmute(drug = A,
              MIC_S = B,
              MIC_R = C,
              disk_dose = ifelse(has_zone_diameters, E, NA_character_),
              disk_S = ifelse(has_zone_diameters, `F`, NA_character_),
              disk_R = ifelse(has_zone_diameters, G, NA_character_)) %>% 
    filter(!is.na(drug),
           !(is.na(MIC_S) & is.na(MIC_R) & is.na(disk_S) & is.na(disk_R)),
           MIC_S %unlike% "(MIC|S ≤|note)",
           MIC_S %unlike% "^[-]",
           drug != MIC_S,) %>% 
    mutate(administration = case_when(drug %like% "[( ]oral" ~ "oral",
                                      drug %like% "[( ]iv"   ~ "iv",
                                      TRUE ~ NA_character_),
           uti = ifelse(drug %like% "(UTI|urinary|urine)", TRUE, FALSE),
           systemic = ifelse(drug %like% "(systemic|septic)", TRUE, FALSE),
           mo = ifelse(drug %like% "([.]|spp)", get_mo(drug), mo_sheet)) %>% 
    # clean disk doses
    mutate(disk_dose = clean_character(disk_dose, remove = "[^0-9.-]")) %>% 
    # clean MIC and disk values
    mutate(MIC_S = gsub(".,.", "", MIC_S), # remove superscript notes with comma, like 0.5^2,3
           MIC_R = gsub(".,.", "", MIC_R),
           disk_S = gsub(".,.", "", disk_S),
           disk_R = gsub(".,.", "", disk_R),
           MIC_S = clean_double(MIC_S), # make them valid numeric values
           MIC_R = clean_double(MIC_R),
           disk_S = clean_integer(disk_S),
           disk_R = clean_integer(disk_R),
           # invalid MIC values have a superscript text, delete those
           MIC_S = ifelse(has_superscript(MIC_S),
                          substr(MIC_S, 1, nchar(MIC_S) - 1),
                          MIC_S),
           MIC_R = ifelse(has_superscript(MIC_R),
                          substr(MIC_R, 1, nchar(MIC_R) - 1),
                          MIC_R),
           # and some are just awful
           MIC_S = ifelse(MIC_S == 43.4, 4, MIC_S),
           MIC_R = ifelse(MIC_R == 43.4, 4, MIC_R),
    ) %>% 
    # clean drug names
    mutate(drug = gsub(" ?[(, ].*$", "", drug),
           drug = gsub("[1-9]+$", "", drug),
           ab = as.ab(drug)) %>% 
    select(ab, mo, everything(), -drug) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  # new row for every different MO mentioned
  for (i in 1:nrow(cleaned)) {
    mo <- cleaned[i, "mo", drop = TRUE]
    if (grepl(pattern = "|", mo, fixed = TRUE)) {
      mo_vect <- unlist(strsplit(mo, "|", fixed = TRUE))
      cleaned[i, "mo"] <- mo_vect[1]
      for (j in seq_len(length(mo_vect))) {
        cleaned <- bind_rows(cleaned, cleaned[i , , drop = FALSE])
        cleaned[nrow(cleaned), "mo"] <- mo_vect[j]
      }
    }
  }
  
  cleaned <- cleaned %>% 
    distinct(ab, mo, administration, uti, systemic, .keep_all = TRUE) %>% 
    arrange(ab, mo) %>% 
    mutate_at(c("MIC_S", "MIC_R", "disk_S", "disk_R"), as.double) %>%
    pivot_longer(c("MIC_S", "MIC_R", "disk_S", "disk_R"), "type") %>% 
    mutate(method = ifelse(type %like% "MIC", "MIC", "DISK"), 
           type = gsub("^.*_", "breakpoint_", type)) %>%
    pivot_wider(names_from = type, values_from = value) %>% 
    mutate(guideline = guideline_name,
           disk_dose = ifelse(method == "DISK", disk_dose, NA_character_),
           mo = ifelse(mo == "", mo_sheet, mo)) %>% 
    filter(!(is.na(breakpoint_S) & is.na(breakpoint_R))) %>% 
    # comply with rsi_translation for now
    transmute(guideline,
              method, 
              site = case_when(uti ~ "UTI",
                               systemic ~ "Systemic",
                               TRUE ~ administration),
              mo, ab, 
              ref_tbl = sheet.bak, 
              disk_dose = ifelse(!is.na(disk_dose), paste0(disk_dose, "ug"), NA_character_),
              breakpoint_S, 
              breakpoint_R) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  if (!is.null(uncertainties)) {
    print(uncertainties %>% distinct(input, mo, .keep_all = TRUE))
  }

  message("Estimated: ", probable_rows, ", gained: ", cleaned %>% count(ab) %>% nrow())
  message("------------------------")
  cleaned
}


# Actual import -----------------------------------------------------------

file <- "data-raw/v_11.0_Breakpoint_Tables.xlsx"
sheets <- readxl::excel_sheets(file)
guideline_name <- "EUCAST 2021"

sheets_to_analyse <- sheets[!sheets %in% c("Content", "Changes", "Notes", "Guidance", "Dosages", "Technical uncertainty", "Topical agents")]

# takes the longest time:
new_EUCAST <- read_EUCAST(sheet = sheets_to_analyse[1],
                          file = file, 
                          guideline_name = guideline_name) 
for (i in 2:length(sheets_to_analyse)) {
  tryCatch(
    new_EUCAST <<- bind_rows(new_EUCAST,
                             read_EUCAST(sheet = sheets_to_analyse[i],
                                         file = file, 
                                         guideline_name = guideline_name))
    , error = function(e) message(e$message))
}

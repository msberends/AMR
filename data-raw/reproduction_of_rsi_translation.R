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

library(dplyr)
library(readr)
library(tidyr)

# Installed WHONET software on Windows (http://www.whonet.org/software.html),
#    imported C:\WHONET\Codes\DRGLST1.txt
DRGLST1 <- readr::read_tsv("data-raw/DRGLST1.txt", na = c("", "NA", "-"))
rsi_trans <- DRGLST1 %>%
  # only keep CLSI and EUCAST guidelines:
  filter(GUIDELINES %like% "^(CLSI|EUCST)")
if (any(is.na(rsi_trans$BREAKPOINT_TYPE)) | !"Human" %in% rsi_trans$BREAKPOINT_TYPE) {
  stop("Check column BREAKPOINT_TYPE - something is WRONG!")
}
sort(unique(rsi_trans$GUIDELINES))
rsi_trans <- rsi_trans %>% 
  ##### If looking for adding a specific guideline, do it here!
  filter(GUIDELINES == "CLSI21") %>% 
  #####
  filter(BREAKPOINT_TYPE == "Human") %>% 
  mutate(DISK_S = ifelse(as.double(DISK_S) > 50, 50, DISK_S),
         MIC_R = ifelse(as.double(MIC_R) %in% c(1025, 129, 513), as.double(MIC_R) - 1, MIC_R)) %>%
  # set a nice layout:
  transmute(guideline = gsub("([0-9]+)$", " 20\\1", gsub("EUCST", "EUCAST", GUIDELINES)),
            method = TESTMETHOD,
            site = SITE_INF,
            mo = as.mo(ORG_CODE),
            ab = as.ab(WHON5_CODE),
            ref_tbl = REF_TABLE,
            dose_disk = POTENCY,
            S_disk = as.disk(DISK_S),
            R_disk = as.disk(DISK_R),
            S_mic = as.mic(MIC_S),
            R_mic = as.mic(MIC_R)) %>%
  filter(!is.na(mo),
         !is.na(ab),
         !mo %in% c("UNKNOWN", "B_GRAMN", "B_GRAMP", "F_FUNGUS", "F_YEAST")) %>%
  arrange(desc(guideline), mo, ab)

print(mo_failures())

# create 2 tables: MIC and disk
tbl_mic <- rsi_trans %>%
  filter(method == "MIC") %>%
  mutate(breakpoint_S = as.double(S_mic), breakpoint_R = as.double(R_mic))
tbl_disk <- rsi_trans %>%
  filter(method == "DISK") %>%
  mutate(breakpoint_S = as.double(S_disk), breakpoint_R = as.double(R_disk))

# merge them so every record is a unique combination of method, mo and ab
rsi_trans <- bind_rows(tbl_mic, tbl_disk) %>%
  rename(disk_dose = dose_disk) %>% 
  mutate(disk_dose = gsub("µ", "u", disk_dose)) %>% 
  select(-ends_with("_mic"), -ends_with("_disk"))

# add extra CLSI general guidelines
# Installed WHONET software on Windows (http://www.whonet.org/software.html),
#    imported C:\WHONET\Codes\DRGLST.txt
clsi_general <- readr::read_tsv("data-raw/DRGLST.txt") %>%
  filter(CLSI == "X") %>%
  select(WHON5_CODE, 
         disk_dose = POTENCY, 
         starts_with("CLSI"), 
         -c(CLSI, CLSI_ORDER)) %>%
  mutate_at(vars(matches("CLSI")), as.double) %>%
  pivot_longer(-c(WHON5_CODE, disk_dose)) %>%
  mutate(method = ifelse(name %like% "_D", "DISK", "MIC"),
         breakpoint = paste0("breakpoint_", gsub(".*([A-Z])$", "\\1", name)), 
         guideline = paste0("CLSI 20", cleaner::clean_integer(name))) %>%
  filter(breakpoint != "breakpoint_I", !is.na(value)) %>%
  select(-name) %>%
  pivot_wider(names_from = breakpoint, values_from = value) %>% 
  transmute(guideline, 
            method, 
            site = NA_character_, 
            mo = as.mo("UNKNOWN"),
            ab = as.ab(WHON5_CODE),
            ref_tbl = "Generic CLSI rules", 
            disk_dose = gsub("/", "-", disk_dose, fixed = TRUE), 
            breakpoint_S, 
            breakpoint_R)


# add new EUCAST with read_EUCAST.R

# 2020-04-14 did that now for 2019 and 2020
rsi_trans <- rsi_trans %>%
  filter(guideline != "EUCAST 2019") %>% 
  bind_rows(new_EUCAST) %>% 
  bind_rows(clsi_general) %>% 
  mutate(uti = site %like% "(UTI|urinary|urine)") %>% 
  as.data.frame(stringsAsFactors = FALSE) %>%
  # force classes again
  mutate(mo = as.mo(mo),
         ab = as.ab(ab)) %>% 
  arrange(desc(guideline), ab, mo, method)

# 2021-01-12 did that now for 2021
rsi_trans <- rsi_trans %>%
  mutate(mo = as.character(mo)) %>% 
  bind_rows(new_EUCAST) %>% 
  mutate(uti = site %like% "(UTI|urinary)") %>% 
  as.data.frame(stringsAsFactors = FALSE) %>%
  # force classes again
  mutate(mo = as.mo(mo),
         ab = as.ab(ab)) %>% 
  arrange(desc(guideline), ab, mo, method)

# save to package
rsi_translation <- rsi_trans
usethis::use_data(rsi_translation, overwrite = TRUE)
rm(rsi_trans)
rm(rsi_translation)
devtools::load_all(".")

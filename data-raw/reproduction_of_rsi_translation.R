library(dplyr)

# Installed WHONET 2019 software on Windows (http://www.whonet.org/software.html),
#    opened C:\WHONET\Codes\WHONETCodes.mdb in MS Access
#    and exported table 'DRGLST1' to MS Excel
DRGLST1 <- readxl::read_excel("data-raw/DRGLST1.xlsx", na = c("", "NA", "-"))
rsi_translation <- DRGLST1 %>%
  # only keep CLSI and EUCAST guidelines:
  filter(GUIDELINES %like% "^(CLSI|EUCST)") %>%
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
tbl_mic <- rsi_translation %>%
  filter(method == "MIC") %>%
  mutate(breakpoint_S = as.double(S_mic), breakpoint_R = as.double(R_mic))
tbl_disk <- rsi_translation %>%
  filter(method == "DISK") %>%
  mutate(breakpoint_S = as.double(S_disk), breakpoint_R = as.double(R_disk))

# merge them so every record is a unique combination of method, mo and ab
rsi_translation <- bind_rows(tbl_mic, tbl_disk) %>%
  rename(disk_dose = dose_disk) %>% 
  mutate(disk_dose = gsub("µ", "u", disk_dose)) %>% 
  select(-ends_with("_mic"), -ends_with("_disk"))

# add extra CLSI general guidelines
clsi_general <- read_tsv("data-raw/DRGLST.txt") %>%
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
rsi_translation <- rsi_translation %>%
  # filter(guideline != "EUCAST 2019") %>% 
  bind_rows(new_EUCAST) %>% 
  bind_rows(clsi_general) %>% 
  mutate(uti = site %like% "(UTI|urinary)") %>% 
  as.data.frame(stringsAsFactors = FALSE) %>%
  # force classes again
  mutate(mo = as.mo(mo),
         ab = as.ab(ab)) %>% 
  arrange(desc(guideline), ab, mo, method)

# save to package
usethis::use_data(rsi_translation, overwrite = TRUE)
rm(rsi_translation)
devtools::load_all(".")

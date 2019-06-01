library(dplyr)
library(readxl)

# Installed WHONET 2019 software on Windows (http://www.whonet.org/software.html),
#    opened C:\WHONET\Codes\WHONETCodes.mdb in MS Access
#    and exported table 'DRGLST1' to MS Excel
DRGLST1 <- read_excel("DRGLST1.xlsx")
rsi_translation <- DRGLST1 %>%
  # only keep CLSI and EUCAST guidelines:
  filter(GUIDELINES %like% "^(CLSI|EUCST)") %>%
  # set a nice layout:
  transmute(guideline = gsub("([0-9]+)$", " 20\\1", gsub("EUCST", "EUCAST", GUIDELINES)),
            method = TESTMETHOD,
            mo = as.mo(ORG_CODE),
            ab = as.ab(WHON5_CODE),
            ref_tbl = REF_TABLE,
            dose_disk = POTENCY,
            S_disk = as.disk(DISK_S),
            R_disk = as.disk(DISK_R),
            S_mic = as.mic(MIC_S),
            R_mic = as.mic(MIC_R)) %>%
  filter(!is.na(mo) & !is.na(ab)) %>%
  arrange(desc(guideline), mo, ab)

# create 2 tables: MIC and disk
tbl_mic <- rsi_translation %>%
  filter(method == "MIC") %>%
  select(-ends_with("_disk")) %>%
  mutate(joinstring = paste(guideline, mo, ab))
tbl_disk <- rsi_translation %>%
  filter(method == "DISK") %>%
  select(-S_mic, -R_mic) %>%
  mutate(joinstring = paste(guideline, mo, ab)) %>%
  select(joinstring, ends_with("_disk"))

# merge them so every record is a unique combination of method, mo and ab
rsi_translation <- tbl_mic %>%
  left_join(tbl_disk,
            by = "joinstring") %>%
  select(-joinstring, -method) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  # force classes again
  mutate(mo = as.mo(mo),
         ab = as.ab(ab),
         S_mic = as.mic(S_mic),
         R_mic = as.mic(R_mic),
         S_disk = as.disk(S_disk),
         R_disk = as.disk(R_disk))

# save to package
usethis::use_data(rsi_translation, overwrite = TRUE)
rm(rsi_translation)

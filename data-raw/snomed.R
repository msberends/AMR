# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

library(AMR)
library(tidyverse)

# go to https://www.nictiz.nl/standaardisatie/terminologiecentrum/referentielijsten/micro-organismen/
# read the table from clipboard
snomed <- clipr::read_clip_tbl()
# snomed <- snomed %>%
#   transmute(fullname = trimws(gsub("^genus", "", Omschrijving, ignore.case = TRUE)),
#             snomed = as.integer(Id))
snomed <- snomed %>%
  transmute(fullname = mo_name(Omschrijving),
            snomed = as.integer(Id)) %>% 
  filter(!fullname %like% "unknown")
snomed_trans <- snomed %>%
  group_by(fullname) %>%
  mutate(snomed_list = list(snomed)) %>%
  ungroup() %>%
  select(fullname, snomed = snomed_list) %>%
  distinct(fullname, .keep_all = TRUE)

microorganisms <- AMR::microorganisms %>% 
  left_join(snomed_trans)
# remove the NULLs, set to NA
microorganisms$snomed <- lapply(microorganisms$snomed, function(x) if (length(x) == 0)  NA else x)

microorganisms <- dataset_UTF8_to_ASCII(microorganisms)

usethis::use_data(microorganisms, overwrite = TRUE)
rm(microorganisms)

# OLD ---------------------------------------------------------------------

baseUrl <- 'https://browser.ihtsdotools.org/snowstorm/snomed-ct'
edition <- 'MAIN'
version <- '2019-07-31'

microorganisms.snomed <- data.frame(conceptid = character(0),
                                    mo = character(0),
                                    stringsAsFactors = FALSE)
microorganisms$snomed <- ""

# for (i in 1:50) {
for (i in 1:1000) {
  
  if (i %% 10 == 0) {
    cat(paste0(i, " - ", cleaner::percentage(i / nrow(microorganisms)), "\n"))
  }
  
  mo_data <- microorganisms %>% 
    filter(mo == microorganisms$mo[i]) %>% 
    as.list()
  
  if (!mo_data$rank %in% c("genus", "species")) {
    next
  }
  
  searchTerm <- paste0(
    ifelse(mo_data$rank == "genus", "Genus ", ""),
    mo_data$fullname, 
    " (organism)")
  
  url <- paste0(baseUrl, '/browser/',
                edition, '/', 
                version, 
                '/descriptions?term=', curl::curl_escape(searchTerm),
                '&mode=fullText&activeFilter=true&limit=', 250)
  results <- url %>% 
    httr::GET() %>%
    httr::content(type = "text", encoding = "UTF-8") %>% 
    jsonlite::fromJSON(flatten = TRUE) %>% 
    .$items
  if (NROW(results) == 0) {
    next
  } else {
    message("Adding ", crayon::italic(mo_data$fullname))
  }
  
  tryCatch(
    microorganisms$snomed[i] <- results %>% filter(term == searchTerm) %>% pull(concept.conceptId),
    error = function(e) invisible()
  )
  
  if (nrow(results) > 1) {
      microorganisms.snomed <- microorganisms.snomed %>% 
        bind_rows(tibble(conceptid = results %>% filter(term != searchTerm) %>% pull(concept.conceptId) %>% unique(),
                         mo = as.character(mo_data$mo)))
  }
}

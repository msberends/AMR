
# WORK IN PROGRESS --------------------------------------------------------

# vector with official names, return vector with CIDs
get_CID <- function(ab) {
  CID <- rep(NA_integer_, length(ab))
  p <- progress_estimated(n = length(ab), min_time = 0)
  for (i in 1:length(ab)) {
    p$tick()$print()

    CID[i] <- tryCatch(
      data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                               ab[i],
                               "/cids/TXT?name_type=complete"),
                        showProgress = FALSE)[[1]][1],
      error = function(e) NA_integer_)
    Sys.sleep(0.2)
  }
  CID
}

# returns vector with synonyms (brand names) for a single CID
get_synonyms <- function(CID, clean = TRUE) {
  p <- progress_estimated(n = length(CID), min_time = 0)
  p$tick()$print()

  synonyms_txt <- tryCatch(
    data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/cid/",
                                           CID,
                                           "/synonyms/TXT"),
                                    sep = "\n",
                                    showProgress = FALSE)[[1]],
    error = function(e) NA_character_)

  if (clean == TRUE) {
    # remove txt between brackets
    synonyms_txt <- trimws(gsub("[(].*[)]", "", gsub("[[].*[]]", "", synonyms_txt)))
    # only length 6 to 20 and no txt with reading marks or numbers
    synonyms_txt <- synonyms_txt[nchar(synonyms_txt) %in% c(6:20)
                                 & !synonyms_txt %like% "[-&{},_0-9]"]
    synonyms_txt <- unlist(strsplit(synonyms_txt,  ";", fixed = TRUE))
  }
  synonyms_txt <- synonyms_txt[tolower(synonyms_txt) %in% unique(tolower(synonyms_txt))]
  sort(synonyms_txt)
}

CIDs <- get_CID(antibiotics$official)
synonyms <- character(length(CIDs))
p <- progress_estimated(n = length(synonyms), min_time = 0)
for (i in 365:length(synonyms)) {
  #p$tick()$print()
  if (!is.na(CIDs[i])) {
    synonyms[i] <- paste(get_synonyms(CIDs[i]), collapse = "|")
  }
}

antibiotics$cid <- CIDs
antibiotics$trade_name <- synonyms

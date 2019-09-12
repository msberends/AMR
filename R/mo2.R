# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

# THIS IS A TEST FUNCTION
as.mo2 <- function(x,
                      Becker = FALSE,
                      Lancefield = FALSE,
                      allow_uncertain = TRUE,
                      reference_df = get_mo_source(),
                      property = "mo",
                      initial_search = TRUE,
                      dyslexia_mode = FALSE,
                      force_mo_history = FALSE,
                      debug = FALSE) {
  
  if (!"AMR" %in% base::.packages()) {
    require("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }
  
  # WHONET: xxx = no growth
  x[tolower(as.character(paste0(x, ""))) %in% c("", "xxx", "na", "nan")] <- NA_character_
  
  if (initial_search == TRUE) {
    options(mo_failures = NULL)
    options(mo_uncertainties = NULL)
    options(mo_renamed = NULL)
  }
  options(mo_renamed_last_run = NULL)
  
  if (NCOL(x) == 2) {
    # support tidyverse selection like: df %>% select(colA, colB)
    # paste these columns together
    x_vector <- vector("character", NROW(x))
    for (i in 1:NROW(x)) {
      x_vector[i] <- paste(pull(x[i,], 1), pull(x[i,], 2), sep = " ")
    }
    x <- x_vector
  } else {
    if (NCOL(x) > 2) {
      stop('`x` can be 2 columns at most', call. = FALSE)
    }
    x[is.null(x)] <- NA
    
    # support tidyverse selection like: df %>% select(colA)
    if (!is.vector(x) & !is.null(dim(x))) {
      x <- pull(x, 1)
    }
  }
  
  notes <- character(0)
  uncertainties <- data.frame(uncertainty = integer(0),
                              input = character(0),
                              fullname = character(0),
                              renamed_to = character(0),
                              mo = character(0), 
                              stringsAsFactors = FALSE)
  failures <- character(0)
  uncertainty_level <- translate_allow_uncertain(allow_uncertain)
  
  
  # x_input <- x
  # already strip leading and trailing spaces
  #x <- trimws(x, which = "both")
  # only check the uniques, which is way faster
  #x <- unique(x)
  # remove empty values (to later fill them in again with NAs)
  # ("xxx" is WHONET code for 'no growth')
  # x <- x[!is.na(x)
  #        & !is.null(x)
  #        & !identical(x, "")
  #        & !identical(x, "xxx")]
  
  # conversion of old MO codes from v0.5.0 (ITIS) to later versions (Catalogue of Life)
  if (any(x %like% "^[BFP]_[A-Z]{3,7}") & !all(x %in% microorganisms$mo)) {
    leftpart <- gsub("^([BFP]_[A-Z]{3,7}).*", "\\1", x)
    if (any(leftpart %in% names(mo_codes_v0.5.0))) {
      rightpart <- gsub("^[BFP]_[A-Z]{3,7}(.*)", "\\1", x)
      leftpart <- mo_codes_v0.5.0[leftpart]
      x[!is.na(leftpart)] <- paste0(leftpart[!is.na(leftpart)], rightpart[!is.na(leftpart)])
    }
    # now check if some are still old
    still_old <- x[x %in% names(mo_codes_v0.5.0)]
    if (length(still_old) > 0) {
      x[x %in% names(mo_codes_v0.5.0)] <- data.frame(old = still_old, stringsAsFactors = FALSE) %>%
        left_join(data.frame(old = names(mo_codes_v0.5.0),
                             new = mo_codes_v0.5.0,
                             stringsAsFactors = FALSE), by = "old") %>%
        # if they couldn't be found, replace them with the old ones again,
        # so they will throw a warning in the end
        mutate(new = ifelse(is.na(new), old, new)) %>%
        pull(new)
    }
  }
  
  # # defined df to check for
  # if (!is.null(reference_df)) {
  #   if (!mo_source_isvalid(reference_df)) {
  #     stop("`reference_df` must contain a column `mo` with values from the 'microorganisms' data set.", call. = FALSE)
  #   }
  #   reference_df <- reference_df %>% filter(!is.na(mo))
  #   # keep only first two columns, second must be mo
  #   if (colnames(reference_df)[1] == "mo") {
  #     reference_df <- reference_df[, c(2, 1)]
  #   } else {
  #     reference_df <- reference_df[, c(1, 2)]
  #   }
  #   colnames(reference_df)[1] <- "x"
  #   # remove factors, just keep characters
  #   suppressWarnings(
  #     reference_df[] <- lapply(reference_df, as.character)
  #   )
  # }
  # 
  # # all empty
  # if (all(identical(trimws(x_input), "") | is.na(x_input) | length(x) == 0)) {
  #   if (property == "mo") {
  #     return(to_class_mo(rep(NA_character_, length(x_input))))
  #   } else {
  #     return(rep(NA_character_, length(x_input)))
  #   }
  #   
  # } else if (all(x %in% reference_df[, 1][[1]])) {
  #   # all in reference df
  #   colnames(reference_df)[1] <- "x"
  #   suppressWarnings(
  #     x <- data.frame(x = x, stringsAsFactors = FALSE) %>%
  #       left_join(reference_df, by = "x") %>%
  #       left_join(AMR::microorganisms, by = "mo") %>%
  #       pull(property)
  #   )
  #   
  # } else if (all(x %in% AMR::microorganisms$mo)) {
  #   # existing mo codes when not looking for property "mo", like mo_genus("B_ESCHR_COL")
  #   y <- microorganismsDT[prevalence == 1][data.table(mo = x), on = "mo", ..property][[1]]
  #   if (any(is.na(y))) {
  #     y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(mo = x[is.na(y)]),
  #                                                      on = "mo",
  #                                                      ..property][[1]]
  #   }
  #   if (any(is.na(y))) {
  #     y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(mo = x[is.na(y)]),
  #                                                      on = "mo",
  #                                                      ..property][[1]]
  #   }
  #   x <- y
  #   
  # # } else if (all(x %in% read_mo_history(uncertainty_level,
  # #                                       force = force_mo_history)$x)) {
  # #   # previously found code
  # #   x <- microorganismsDT[data.table(mo = get_mo_history(x,
  # #                                                        uncertainty_level,
  # #                                                        force = force_mo_history)),
  # #                         on = "mo", ..property][[1]]
  #   
  # } else if (all(tolower(x) %in% microorganismsDT$fullname_lower)) {
  #   # we need special treatment for very prevalent full names, they are likely!
  #   # e.g. as.mo("Staphylococcus aureus")
  #   y <- microorganismsDT[prevalence == 1][data.table(fullname_lower = tolower(x)), on = "fullname_lower", ..property][[1]]
  #   if (any(is.na(y))) {
  #     y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(fullname_lower = tolower(x[is.na(y)])),
  #                                                      on = "fullname_lower",
  #                                                      ..property][[1]]
  #   }
  #   if (any(is.na(y))) {
  #     y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(fullname_lower = tolower(x[is.na(y)])),
  #                                                      on = "fullname_lower",
  #                                                      ..property][[1]]
  #   }
  #   x <- y
  #   
  # } else if (all(toupper(x) %in% AMR::microorganisms.codes$code)) {
  #   # commonly used MO codes
  #   y <- as.data.table(AMR::microorganisms.codes)[data.table(code = toupper(x)), on = "code", ]
  #   # save them to history
  #   set_mo_history(x, y$mo, 0, force = force_mo_history)
  #   
  #   x <- microorganismsDT[data.table(mo = y[["mo"]]), on = "mo", ..property][[1]]
  #   
  # } else if (!all(x %in% AMR::microorganisms[, property])) {
  #   
  if (1 == 1) {
    strip_whitespace <- function(x) {
      # all whitespaces (tab, new lines, etc.) should be one space
      # and spaces before and after should be omitted
      trimws(gsub("[\\s]+", " ", x, perl = TRUE), which = "both")
    }
    
    
    x_new <- rep(NA_character_, length(x))
    
    # keep only dots, letters, numbers, slashes, spaces and dashes
    x <- gsub("[^.a-zA-Z0-9/ \\-]+", "", x)
    # remove spp and species
    x <- gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x, ignore.case = TRUE)
    # remove 'genus' as first word
    x <- gsub("^genus ", "", x, ignore.case = TRUE)
    # remove 'uncertain'-like texts
    x <- trimws(gsub("(uncertain|susp[ie]c[a-z]+|verdacht)", "", x, ignore.case = TRUE))
    x <- strip_whitespace(x)
    x_backup <- x
    
    # remove spp and species
    #x <- gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x_backup, ignore.case = TRUE)
    #x <- strip_whitespace(x)
    
    x_backup_without_spp <- x
    x_species <- paste(x, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(gruppe|groep|grupo|gruppo|groupe)", "group", x, ignore.case = TRUE)
    x <- gsub("(vergroen)[a-z]*", "viridans", x, ignore.case = TRUE)
    x <- gsub("(hefe|gist|gisten|levadura|lievito|fermento|levure)[a-z]*", "yeast", x, ignore.case = TRUE)
    x <- gsub("(schimmels?|mofo|molde|stampo|moisissure|fungi)[a-z]*", "fungus", x, ignore.case = TRUE)
    x <- gsub("Fungus[ph|f]rya", "Fungiphrya", x, ignore.case = TRUE)
    # remove non-text in case of "E. coli" except dots and spaces
    #   x <- gsub("[^.a-zA-Z0-9/ \\-]+", "", x)
    # replace minus by a space
    x <- gsub("-+", " ", x)
    # replace hemolytic by haemolytic
    x <- gsub("ha?emoly", "haemoly", x, ignore.case = TRUE)
    # place minus back in streptococci
    x <- gsub("(alpha|beta|gamma).?ha?emoly", "\\1-haemoly", x, ignore.case = TRUE)
    # remove genus as first word
    #  x <- gsub("^genus ", "", x, ignore.case = TRUE)
    # remove 'uncertain' like texts
    # x <- trimws(gsub("(uncertain|susp[ie]c[a-z]+|verdacht)", "", x, ignore.case = TRUE))
    # allow characters that resemble others = dyslexia_mode ----
    if (dyslexia_mode == TRUE) {
      x <- tolower(x)
      x <- gsub("[iy]+", "[iy]+", x)
      x <- gsub("(c|k|q|qu|s|z|x|ks)+", "(c|k|q|qu|s|z|x|ks)+", x)
      x <- gsub("(ph|hp|f|v)+", "(ph|hp|f|v)+", x)
      x <- gsub("(th|ht|t)+", "(th|ht|t)+", x)
      x <- gsub("a+", "a+", x)
      x <- gsub("u+", "u+", x)
      # allow any ending of -um, -us, -ium, -icum, -ius, -icus, -ica and -a (needs perl for the negative backward lookup):
      x <- gsub("(u\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+a\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+u\\+m)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("e+", "e+", x, ignore.case = TRUE)
      x <- gsub("o+", "o+", x, ignore.case = TRUE)
      x <- gsub("(.)\\1+", "\\1+", x)
      # allow ending in -en or -us
      x <- gsub("e\\+n(?![a-z[])", "(e+n|u+(c|k|q|qu|s|z|x|ks)+)", x, ignore.case = TRUE, perl = TRUE)
      # if the input is longer than 10 characters, allow any constant between all characters, as some might have forgotten a character
      # this will allow "Pasteurella damatis" to be correctly read as "Pasteurella dagmatis".
      constants <- paste(letters[!letters %in% c("a", "e", "i", "o", "u")], collapse = "")
      #x[nchar(x_backup_without_spp) > 10] <- gsub("([a-z])([a-z])", paste0("\\1[", constants, "]?\\2"), x[nchar(x_backup_without_spp) > 10], ignore.case = TRUE)
      x[nchar(x_backup_without_spp) > 10] <- gsub("[+]", paste0("+[", constants, "]?"), x[nchar(x_backup_without_spp) > 10])
    }
    x <- strip_whitespace(x)
    
    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    x_trimmed_without_group <- gsub(" gro.u.p$", "", x_trimmed, ignore.case = TRUE)
    # remove last part from "-" or "/"
    x_trimmed_without_group <- gsub("(.*)[-/].*", "\\1", x_trimmed_without_group)
    # replace space and dot by regex sign
    x_withspaces <- gsub("[ .]+", ".* ", x)
    x <- gsub("[ .]+", ".*", x)
    # add start en stop regex
    x <- paste0('^', x, '$')
    x_withspaces_start_only <- paste0('^', x_withspaces)
    x_withspaces_end_only <- paste0(x_withspaces, '$')
    x_withspaces_start_end <- paste0('^', x_withspaces, '$')
    
    if (isTRUE(debug)) {
      print(data.frame(
        x_backup,
        x,
        x_species,
        x_withspaces_start_only,
        x_withspaces_end_only,
        x_withspaces_start_end,
        x_backup_without_spp,
        x_trimmed,
        x_trimmed_species,
        x_trimmed_without_group), right = FALSE)
      # cat(paste0('x                       "', x, '"\n'))
      # cat(paste0('x_species               "', x_species, '"\n'))
      # cat(paste0('x_withspaces_start_only "', x_withspaces_start_only, '"\n'))
      # cat(paste0('x_withspaces_end_only   "', x_withspaces_end_only, '"\n'))
      # cat(paste0('x_withspaces_start_end  "', x_withspaces_start_end, '"\n'))
      # cat(paste0('x_backup                "', x_backup, '"\n'))
      # cat(paste0('x_backup_without_spp    "', x_backup_without_spp, '"\n'))
      # cat(paste0('x_trimmed               "', x_trimmed, '"\n'))
      # cat(paste0('x_trimmed_species       "', x_trimmed_species, '"\n'))
      # cat(paste0('x_trimmed_without_group "', x_trimmed_without_group, '"\n'))
    }
    
    #progress <- progress_estimated(n = length(x), min_time = 3)
    
    # THE NEW WAY ----
    nothing_more_to_do <- function() !any(is.na(x_new) & !is.na(x_backup))
    
    lookup_regexes <- function(data, property, regex) {
      prop <- regex %>% 
        sapply(function(pattern) pull(data, property) %like% pattern) %>%
        as.data.frame() %>% 
        lapply(function(c) suppressWarnings(min(pull(data, property)[c]))) %>% 
        unlist()
      if (is.null(prop)) {
        return(rep(NA, length(regex)))
      }
      DT <- data.table(prop)
      colnames(DT) <- property
      microorganismsDT[DT, on = property, "mo"][[1]]
    }
    
    # LATER: only unique if more than 500 values, of which max 85% distinct?
    
    x_backup_upper <- toupper(x_backup)
    x_backup_lower <- tolower(x_backup)
    
    # exclude all viruses (there is no fullname containing 'virus' in the data set)
    x_new[x_backup %like% "virus"] <- "UNKNOWN"
    
    # try available fields in the microorganisms data set
    x_new[x_backup_upper %in% microorganisms$mo] <- microorganismsDT[data.table(mo = x_backup_upper[x_backup_upper %in% microorganisms$mo]), on = "mo", "mo"][[1]]
    x_new[x_backup_lower %in% microorganismsDT$fullname_lower] <- microorganismsDT[data.table(fullname_lower = x_backup_lower[x_backup_lower %in% microorganismsDT$fullname_lower]), on = "fullname_lower", "mo"][[1]]
    # x_new[x_backup %in% microorganisms$col_id] <- microorganismsDT[data.table(col_id = as.integer(x_backup[x_backup %in% microorganisms$col_id])), on = "col_id", ..property][[1]]
    
    # old names
    old_names <- x_backup[x_backup_lower %in% microorganisms.oldDT$fullname_lower]
    x_new[x_backup_lower %in% microorganisms.oldDT$fullname_lower] <- microorganismsDT[microorganisms.oldDT[data.table(fullname_lower = x_backup_lower[x_backup_lower %in% microorganisms.oldDT$fullname_lower]), on = "fullname_lower", "col_id_new"], on = c("col_id" = "col_id_new"), "mo"][[1]]
    
    if (nothing_more_to_do()) {
      if (property != "mo") {
        return(microorganismsDT[data.table(mo = x_new), on = "mo", ..property][[1]])
      } else {
        return(to_class_mo(x_new))
      }
    }
    
    # codes from the microorganisms.codes data set
    x_new[x_backup_upper %in% microorganisms.codes$code] <- as.data.table(microorganisms.codes)[data.table(code = x_backup_upper[x_backup_upper %in% microorganisms.codes$code]), on = "code", "mo"][[1]]
    if (!is.null(reference_df)) {
      colnames(reference_df)[1] <- "code"
      x_new[x_backup_upper %in% reference_df$code] <- as.data.table(reference_df)[data.table(code = x_backup_upper[x_backup_upper %in% reference_df$code]), on = "code", mo][[1]]
      if (!all(x_new %in% microorganisms$mo, na.rm = TRUE)) {
        warning("Values ", paste(x_new[!x_new %in% c(NA, microorganisms$mo)], collapse = ", "), " found in reference_df, but these are not valid MO codes.", call. = FALSE)
        x_new[!x_new %in% c(NA, microorganisms$mo)] <- "UNKNOWN"
      }
    }
    
    x_new[x_backup_upper %in% c("MRSA", "MSSA", "VISA", "VRSA")] <- "B_STPHY_AUR"
    x_new[x_backup_upper %in% c("MRSE", "MSSE")] <- "B_STPHY_EPI"
    x_new[x_backup_upper %in% c("AIEC", "ATEC", "DAEC", "EAEC", "EHEC", "EIEC", "EPEC", "ETEC", "NMEC", "STEC", "UPEC")
          | x_backup_upper %like% "O?(26|103|104|111|121|145|157)"] <- "B_ESCHR_COL"
    x_new[x_backup_upper %in% c("MRPA")] <- "B_PSDMN_AER"
    x_new[x_backup_upper %in% c("CRSM")] <- "B_STNTR_MAL"
    x_new[x_backup_upper %in% c("PISP", "PRSP", "VISP", "VRSP")] <- "B_STRPT_PNE"
    x_new[x_backup_upper %in% c("VRE")
          | x_backup %like% "(enterococci|enterokok|enterococo)[a-z]*?$"] <- "B_ENTRC"
    
    # start showing progress bar here
    progress <- progress_estimated(n = 3, min_time = 0)
    # most prevalent (1)
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", x_withspaces_start_end[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", paste0(trimws(x_withspaces_start_only), " ")[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", paste0(" ", trimws(x_withspaces_end_only))[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", x_trimmed[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", x_trimmed_without_group[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 1] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new)])
    if (nothing_more_to_do()) {
      if (property != "mo") {
        return(microorganismsDT[data.table(mo = x_new), on = "mo", ..property][[1]])
      } else {
        return(to_class_mo(x_new))
      }
    }
    progress$tick()$print()
    # less prevalent (2)
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", x_withspaces_start_end[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", paste0(trimws(x_withspaces_start_only), " ")[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", paste0(" ", trimws(x_withspaces_end_only))[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", x_trimmed[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", x_trimmed_without_group[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 2] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new)])
    if (nothing_more_to_do()) {
      if (property != "mo") {
        return(microorganismsDT[data.table(mo = x_new), on = "mo", ..property][[1]])
      } else {
        return(to_class_mo(x_new))
      }
    }
    progress$tick()$print()
    # least prevalent (3)
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", x_withspaces_start_end[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", paste0(trimws(x_withspaces_start_only), " ")[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", paste0(" ", trimws(x_withspaces_end_only))[!is.na(x_backup) & is.na(x_new)])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", x_trimmed[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", x_trimmed_without_group[!is.na(x_backup) & is.na(x_new) & nchar(x_backup) >= 6])
    x_new[!is.na(x_backup) & is.na(x_new)] <- microorganismsDT[prevalence == 3] %>% lookup_regexes("fullname", x_withspaces_start_only[!is.na(x_backup) & is.na(x_new)])
    
    # all others are UNKNOWN
    x_new[!is.na(x_backup) & is.na(x_new)] <- "UNKNOWN"
    progress$tick()$print()
    
    return(to_class_mo(x_new))
    #for (i in 1:length(x)) {
    for (i in character(0)) {
      
      x[i] <- "UNKNOWN"
      next
      
      # progress$tick()$print()
      
      # if (initial_search == TRUE) {
      #   found <- microorganismsDT[mo == get_mo_history(x_backup[i],
      #                                                  uncertainty_level,
      #                                                  force = force_mo_history),
      #                             ..property][[1]]
      #   # previously found result
      #   if (length(found) > 0) {
      #     x[i] <- found[1L]
      #     next
      #   }
      # }
      
      found <- microorganismsDT[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid MO code
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      
      found <- microorganismsDT[fullname_lower %in% tolower(c(x_backup[i], x_backup_without_spp[i])), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      found <- microorganismsDT[col_id == x_backup[i], ..property][[1]]
      # is a valid Catalogue of Life ID
      if (NROW(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      
      # WHONET: xxx = no growth
      if (tolower(as.character(paste0(x_backup_without_spp[i], ""))) %in% c("", "xxx", "na", "nan")) {
        x[i] <- NA_character_
        next
      }
      
      if (tolower(x_backup_without_spp[i]) %in% c("other", "none", "unknown")) {
        # empty and nonsense values, ignore without warning
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      # check for very small input, but ignore the O antigens of E. coli
      if (nchar(gsub("[^a-zA-Z]", "", x_trimmed[i])) < 3
          & !x_backup_without_spp[i] %like% "O?(26|103|104|104|111|121|145|157)") {
        # check if search term was like "A. species", then return first genus found with ^A
        # if (x_backup[i] %like% "[a-z]+ species" | x_backup[i] %like% "[a-z] spp[.]?") {
        #   # get mo code of first hit
        #   found <- microorganismsDT[fullname %like% x_withspaces_start_only[i], mo]
        #   if (length(found) > 0) {
        #     mo_code <- found[1L] %>% strsplit("_") %>% unlist() %>% .[1:2] %>% paste(collapse = "_")
        #     found <- microorganismsDT[mo == mo_code, ..property][[1]]
        #     # return first genus that begins with x_trimmed, e.g. when "E. spp."
        #     if (length(found) > 0) {
        #       x[i] <- found[1L]
        #       if (initial_search == TRUE) {
        #         set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        #       }
        #       next
        #     }
        #   }
        # }
        # fewer than 3 chars and not looked for species, add as failure
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      if (x_backup_without_spp[i] %like% "virus") {
        # there is no fullname like virus, so don't try to coerce it
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      # translate known trivial abbreviations to genus + species ----
      if (!is.na(x_trimmed[i])) {
        if (toupper(x_backup_without_spp[i]) %in% c('MRSA', 'MSSA', 'VISA', 'VRSA')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_AUR', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c('MRSE', 'MSSE')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_EPI', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == "VRE"
            | x_backup_without_spp[i] %like% '(enterococci|enterokok|enterococo)[a-z]*?$')  {
          x[i] <- microorganismsDT[mo == 'B_ENTRC', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        # support for:
        # - AIEC (Adherent-Invasive E. coli)
        # - ATEC (Atypical Entero-pathogenic E. coli)
        # - DAEC (Diffusely Adhering E. coli)
        # - EAEC (Entero-Aggresive E. coli)
        # - EHEC (Entero-Haemorrhagic E. coli)
        # - EIEC (Entero-Invasive E. coli)
        # - EPEC (Entero-Pathogenic E. coli)
        # - ETEC (Entero-Toxigenic E. coli)
        # - NMEC (Neonatal Meningitis‐causing E. coli)
        # - STEC (Shiga-toxin producing E. coli)
        # - UPEC (Uropathogenic E. coli)
        if (toupper(x_backup_without_spp[i]) %in% c("AIEC", "ATEC", "DAEC", "EAEC", "EHEC", "EIEC", "EPEC", "ETEC", "NMEC", "STEC", "UPEC")
            # also support O-antigens of E. coli: O26, O103, O104, O111, O121, O145, O157
            | x_backup_without_spp[i] %like% "O?(26|103|104|111|121|145|157)") {
          x[i] <- microorganismsDT[mo == 'B_ESCHR_COL', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == 'MRPA') {
          # multi resistant P. aeruginosa
          x[i] <- microorganismsDT[mo == 'B_PSDMN_AER', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == 'CRS'
            | toupper(x_backup_without_spp[i]) == 'CRSM') {
          # co-trim resistant S. maltophilia
          x[i] <- microorganismsDT[mo == 'B_STNTR_MAL', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c('PISP', 'PRSP', 'VISP', 'VRSP')) {
          # peni I, peni R, vanco I, vanco R: S. pneumoniae
          x[i] <- microorganismsDT[mo == 'B_STRPT_PNE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '^G[ABCDFGHK]S$') {
          # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRB)
          x[i] <- microorganismsDT[mo == gsub("G([ABCDFGHK])S", "B_STRPT_GR\\1", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '(streptococ|streptokok).* [ABCDFGHK]$') {
          # Streptococci in different languages, like "estreptococos grupo B"
          x[i] <- microorganismsDT[mo == gsub(".*(streptococ|streptokok|estreptococ).* ([ABCDFGHK])$", "B_STRPT_GR\\2", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'group [ABCDFGHK] (streptococ|streptokok|estreptococ)') {
          # Streptococci in different languages, like "Group A Streptococci"
          x[i] <- microorganismsDT[mo == gsub(".*group ([ABCDFGHK]) (streptococ|streptokok|estreptococ).*", "B_STRPT_GR\\1", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'haemoly.*strept') {
          # Haemolytic streptococci in different languages
          x[i] <- microorganismsDT[mo == 'B_STRPT_HAE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese) ----
        if (x_backup_without_spp[i] %like% '[ck]oagulas[ea] negatie?[vf]'
            | x_trimmed[i] %like% '[ck]oagulas[ea] negatie?[vf]'
            | x_backup_without_spp[i] %like% '[ck]o?ns[^a-z]?$') {
          # coerce S. coagulase negative
          x[i] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '[ck]oagulas[ea] positie?[vf]'
            | x_trimmed[i] %like% '[ck]oagulas[ea] positie?[vf]'
            | x_backup_without_spp[i] %like% '[ck]o?ps[^a-z]?$') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        # streptococcal groups: milleri and viridans
        if (x_trimmed[i] %like% 'strepto.* milleri'
            | x_backup_without_spp[i] %like% 'strepto.* milleri'
            | x_backup_without_spp[i] %like% 'mgs[^a-z]?$') {
          # Milleri Group Streptococcus (MGS)
          x[i] <- microorganismsDT[mo == 'B_STRPT_MIL', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_trimmed[i] %like% 'strepto.* viridans'
            | x_backup_without_spp[i] %like% 'strepto.* viridans'
            | x_backup_without_spp[i] %like% 'vgs[^a-z]?$') {
          # Viridans Group Streptococcus (VGS)
          x[i] <- microorganismsDT[mo == 'B_STRPT_VIR', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'gram[ -]?neg.*'
            | x_backup_without_spp[i] %like% 'negatie?[vf]'
            | x_trimmed[i] %like% 'gram[ -]?neg.*') {
          # coerce Gram negatives
          x[i] <- microorganismsDT[mo == 'B_GRAMN', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'gram[ -]?pos.*'
            | x_backup_without_spp[i] %like% 'positie?[vf]'
            | x_trimmed[i] %like% 'gram[ -]?pos.*') {
          # coerce Gram positives
          x[i] <- microorganismsDT[mo == 'B_GRAMP', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% "salmonella [a-z]+ ?.*") {
          if (x_backup_without_spp[i] %like% "Salmonella group") {
            # Salmonella Group A to Z, just return S. species for now
            x[i] <- microorganismsDT[mo == 'B_SLMNL', ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
          } else if (grepl("[sS]almonella [A-Z][a-z]+ ?.*", x_backup_without_spp[i], ignore.case = FALSE)) {
            # Salmonella with capital letter species like "Salmonella Goettingen" - they're all S. enterica
            x[i] <- microorganismsDT[mo == 'B_SLMNL_ENT', ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
            uncertainties <- rbind(uncertainties,
                                   format_uncertainty_as_df(uncertainty_level = 1,
                                                            input = x_backup_without_spp[i],
                                                            result_mo = "B_SLMNL_ENT"))
          }
          next
        }
        
        # trivial names known to the field:
        if ("meningococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_NESSR_MEN', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if ("gonococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_NESSR_GON', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if ("pneumococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_STRPT_PNE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
      }
      
      # FIRST TRY FULLNAMES AND CODES ----
      # if only genus is available, return only genus
      if (all(!c(x[i], x_trimmed[i]) %like% " ")) {
        found <- microorganismsDT[fullname_lower %in% tolower(c(x_species[i], x_trimmed_species[i])), ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (nchar(x_backup_without_spp[i]) >= 6) {
          found <- microorganismsDT[fullname_lower %like% paste0("^", unregex(x_backup_without_spp[i]), "[a-z]+"), ..property][[1]]
          if (length(found) > 0) {
            x[i] <- found[1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
            next
          }
        }
        # rest of genus only is in allow_uncertain part.
      }
      
      # TRY OTHER SOURCES ----
      # WHONET and other common LIS codes
      if (toupper(x_backup[i]) %in% AMR::microorganisms.codes[, 1]) {
        mo_found <- AMR::microorganisms.codes[toupper(x_backup[i]) == AMR::microorganisms.codes[, 1], "mo"][1L]
        if (length(mo_found) > 0) {
          x[i] <- microorganismsDT[mo == mo_found, ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
      }
      if (!is.null(reference_df)) {
        # self-defined reference
        if (x_backup[i] %in% reference_df[, 1]) {
          ref_mo <- reference_df[reference_df[, 1] == x_backup[i], "mo"]
          if (ref_mo %in% microorganismsDT[, mo]) {
            x[i] <- microorganismsDT[mo == ref_mo, ..property][[1]][1L]
            next
          } else {
            warning("Value '", x_backup[i], "' was found in reference_df, but '", ref_mo, "' is not a valid MO code.", call. = FALSE)
          }
        }
      }
      
      # allow no codes less than 4 characters long, was already checked for WHONET above
      if (nchar(x_backup_without_spp[i]) < 4) {
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      check_per_prevalence <- function(data_to_check,
                                       a.x_backup,
                                       b.x_trimmed,
                                       c.x_trimmed_without_group,
                                       d.x_withspaces_start_end,
                                       e.x_withspaces_start_only,
                                       f.x_withspaces_end_only,
                                       g.x_backup_without_spp) {
        
        # try probable: trimmed version of fullname ----
        found <- data_to_check[fullname_lower %in% tolower(g.x_backup_without_spp), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # try any match keeping spaces ----
        found <- data_to_check[fullname %like% d.x_withspaces_start_end, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not ending with $ ----
        found <- data_to_check[fullname %like% paste0(trimws(e.x_withspaces_start_only), " "), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        found <- data_to_check[fullname %like% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not start with ^ ----
        found <- data_to_check[fullname %like% paste0(" ", trimws(f.x_withspaces_end_only)), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # try a trimmed version
        found <- data_to_check[fullname_lower %like% b.x_trimmed
                               | fullname_lower %like% c.x_trimmed_without_group, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        
        # try splitting of characters in the middle and then find ID ----
        # only when text length is 6 or lower
        # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
        if (nchar(g.x_backup_without_spp) <= 6) {
          x_length <- nchar(g.x_backup_without_spp)
          x_split <- paste0("^",
                            g.x_backup_without_spp %>% substr(1, x_length / 2),
                            '.* ',
                            g.x_backup_without_spp %>% substr((x_length / 2) + 1, x_length))
          found <- data_to_check[fullname %like% x_split, ..property][[1]]
          if (length(found) > 0) {
            return(found[1L])
          }
        }
        
        # try fullname without start and without nchar limit of >= 6 ----
        # like "K. pneu rhino" >> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
        found <- data_to_check[fullname %like% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # didn't found any
        return(NA_character_)
      }
      
      # FIRST TRY VERY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 1],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      # THEN TRY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 2],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      # THEN UNPREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 3],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      # MISCELLANEOUS ----
      
      # look for old taxonomic names ----
      found <- microorganisms.oldDT[fullname_lower == tolower(x_backup[i])
                                    | fullname %like% x_withspaces_start_end[i],]
      if (NROW(found) > 0) {
        col_id_new <- found[1, col_id_new]
        # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
        # mo_ref("Chlamydia psittaci") = "Page, 1968" (with warning)
        # mo_ref("Chlamydophila psittaci") = "Everett et al., 1999"
        if (property == "ref") {
          x[i] <- found[1, ref]
        } else {
          x[i] <- microorganismsDT[col_id == found[1, col_id_new], ..property][[1]]
        }
        options(mo_renamed_last_run = found[1, fullname])
        was_renamed(name_old = found[1, fullname],
                    name_new = microorganismsDT[col_id == found[1, col_id_new], fullname],
                    ref_old = found[1, ref],
                    ref_new = microorganismsDT[col_id == found[1, col_id_new], ref],
                    mo = microorganismsDT[col_id == found[1, col_id_new], mo])
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      
      # check for uncertain results ----
      uncertain_fn <- function(a.x_backup,
                               b.x_trimmed,
                               c.x_withspaces_start_end,
                               d.x_withspaces_start_only,
                               f.x_withspaces_end_only,
                               g.x_backup_without_spp) {
        
        if (uncertainty_level == 0) {
          # do not allow uncertainties
          return(NA_character_)
        }
        
        # UNCERTAINTY LEVEL 1 ----
        if (uncertainty_level >= 1) {
          now_checks_for_uncertainty_level <- 1
          
          # (1) look again for old taxonomic names, now for G. species ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (1) look again for old taxonomic names, now for G. species\n")
          }
          if (isTRUE(debug)) {
            message("Running '", c.x_withspaces_start_end, "' and '", d.x_withspaces_start_only, "'")
          }
          found <- microorganisms.oldDT[fullname %like% c.x_withspaces_start_end
                                        | fullname %like% d.x_withspaces_start_only]
          if (NROW(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
            if (property == "ref") {
              # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
              # mo_ref("Chlamydia psittaci) = "Page, 1968" (with warning)
              # mo_ref("Chlamydophila psittaci) = "Everett et al., 1999"
              x <- found[1, ref]
            } else {
              x <- microorganismsDT[col_id == found[1, col_id_new], ..property][[1]]
            }
            was_renamed(name_old = found[1, fullname],
                        name_new = microorganismsDT[col_id == found[1, col_id_new], fullname],
                        ref_old = found[1, ref],
                        ref_new = microorganismsDT[col_id == found[1, col_id_new], ref],
                        mo = microorganismsDT[col_id == found[1, col_id_new], mo])
            options(mo_renamed_last_run = found[1, fullname])
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = microorganismsDT[col_id == found[1, col_id_new], mo]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(x, property), 1, force = force_mo_history)
            }
            return(x)
          }
          
          # (2) Try with misspelled input ----
          # just rerun with dyslexia_mode = TRUE will used the extensive regex part above
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (2) Try with misspelled input\n")
          }
          if (isTRUE(debug)) {
            message("Running '", a.x_backup, "'")
          }
          # first try without dyslexia mode
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
          if (empty_result(found)) {
            # then with dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
          }
          if (!empty_result(found)) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 1, force = force_mo_history)
            }
            return(found[1L])
          }
        }
        
        # UNCERTAINTY LEVEL 2 ----
        if (uncertainty_level >= 2) {
          now_checks_for_uncertainty_level <- 2
          
          # (3) look for genus only, part of name ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (3) look for genus only, part of name\n")
          }
          if (nchar(g.x_backup_without_spp) > 4 & !b.x_trimmed %like% " ") {
            if (!grepl("^[A-Z][a-z]+", b.x_trimmed, ignore.case = FALSE)) {
              if (isTRUE(debug)) {
                message("Running '", paste(b.x_trimmed, "species"), "'")
              }
              # not when input is like Genustext, because then Neospora would lead to Actinokineospora
              found <- microorganismsDT[fullname_lower %like% paste(b.x_trimmed, "species"), ..property][[1]]
              if (length(found) > 0) {
                x[i] <- found[1L]
                uncertainties <<- rbind(uncertainties,
                                        format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                 input = a.x_backup,
                                                                 result_mo = found_result[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(x, property), 2, force = force_mo_history)
                }
                return(x)
              }
            }
          }
          
          # (4) strip values between brackets ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (4) strip values between brackets\n")
          }
          a.x_backup_stripped <- gsub("( *[(].*[)] *)", " ", a.x_backup)
          a.x_backup_stripped <- trimws(gsub(" +", " ", a.x_backup_stripped))
          if (isTRUE(debug)) {
            message("Running '", a.x_backup_stripped, "'")
          }
          # first try without dyslexia mode
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
          if (empty_result(found)) {
            # then with dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
          }
          if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          
          # (5) inverse input ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (5) inverse input\n")
          }
          a.x_backup_inversed <- paste(rev(unlist(strsplit(a.x_backup, split = " "))), collapse = " ")
          if (isTRUE(debug)) {
            message("Running '", a.x_backup_inversed, "'")
          }
          # first try without dyslexia mode
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
          if (empty_result(found)) {
            # then with dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
          }
          if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          
          # (6) try to strip off half an element from end and check the remains ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (6) try to strip off half an element from end and check the remains\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1) {
            for (i in 1:(length(x_strip) - 1)) {
              lastword <- x_strip[length(x_strip) - i + 1]
              lastword_half <- substr(lastword, 1, as.integer(nchar(lastword) / 2))
              # remove last half of the second term
              x_strip_collapsed <- paste(c(x_strip[1:(length(x_strip) - i)], lastword_half), collapse = " ")
              if (nchar(x_strip_collapsed) >= 4 & nchar(lastword_half) > 2) {
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- microorganismsDT[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
          # (7) try to strip off one element from end and check the remains ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (7) try to strip off one element from end and check the remains\n")
          }
          if (length(x_strip) > 1) {
            for (i in 1:(length(x_strip) - 1)) {
              x_strip_collapsed <- paste(x_strip[1:(length(x_strip) - i)], collapse = " ")
              if (nchar(x_strip_collapsed) >= 6) {
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- microorganismsDT[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
          # (8) check for unknown yeasts/fungi ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (8) check for unknown yeasts/fungi\n")
          }
          if (b.x_trimmed %like% "yeast") {
            found <- "F_YEAST"
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          if (b.x_trimmed %like% "(fungus|fungi)" & !b.x_trimmed %like% "Fungiphrya") {
            found <- "F_FUNGUS"
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                             input = a.x_backup,
                                                             result_mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          # (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome) ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome)\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
            for (i in 2:(length(x_strip))) {
              x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
              if (isTRUE(debug)) {
                message("Running '", x_strip_collapsed, "'")
              }
              # first try without dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
              if (empty_result(found)) {
                # then with dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
              }
              if (!empty_result(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found_result[1L], ..property][[1]]
                # uncertainty level 2 only if searched part contains a space (otherwise it will be found with lvl 3)
                if (x_strip_collapsed %like% " ") {
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
        }
        
        # UNCERTAINTY LEVEL 3 ----
        if (uncertainty_level >= 3) {
          now_checks_for_uncertainty_level <- 3
          
          # (10) try to strip off one element from start and check the remains (any text size) ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (10) try to strip off one element from start and check the remains (any text size)\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
            for (i in 2:(length(x_strip))) {
              x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
              if (isTRUE(debug)) {
                message("Running '", x_strip_collapsed, "'")
              }
              # first try without dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
              if (empty_result(found)) {
                # then with dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
              }
              if (!empty_result(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found, ..property][[1]]
                uncertainties <<- rbind(uncertainties,
                                        format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                 input = a.x_backup,
                                                                 result_mo = found_result[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history)
                }
                return(found[1L])
              }
            }
          }
          # (11) try to strip off one element from end and check the remains (any text size) ----
          # (this is in fact 7 but without nchar limit of >=6)
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (11) try to strip off one element from end and check the remains (any text size)\n")
          }
          if (length(x_strip) > 1) {
            for (i in 1:(length(x_strip) - 1)) {
              x_strip_collapsed <- paste(x_strip[1:(length(x_strip) - i)], collapse = " ")
              if (isTRUE(debug)) {
                message("Running '", x_strip_collapsed, "'")
              }
              # first try without dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug)))
              if (empty_result(found)) {
                # then with dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug)))
              }
              if (!empty_result(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found, ..property][[1]]
                uncertainties <<- rbind(uncertainties,
                                        format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                 input = a.x_backup,
                                                                 result_mo = found_result[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                }
                return(found[1L])
              }
            }
          }
          
          # (12) part of a name (very unlikely match) ----
          if (isTRUE(debug)) {
            cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (12) part of a name (very unlikely match)\n")
          }
          if (isTRUE(debug)) {
            message("Running '", f.x_withspaces_end_only, "'")
          }
          found <- microorganismsDT[fullname %like% f.x_withspaces_end_only]
          if (nrow(found) > 0) {
            found_result <- found[["mo"]]
            if (!empty_result(found_result) & nchar(g.x_backup_without_spp) >= 6) {
              found <- microorganismsDT[mo == found_result[1L], ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history)
              }
              return(found[1L])
            }
          }
        }
        
        # didn't found in uncertain results too
        return(NA_character_)
      }
      x[i] <- uncertain_fn(x_backup[i],
                           x_trimmed[i],
                           x_withspaces_start_end[i],
                           x_withspaces_start_only[i],
                           x_withspaces_end_only[i],
                           x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        # no set_mo_history here - it is already set in uncertain_fn()
        next
      }
      
      # no results found: make them UNKNOWN ----
      x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
      if (initial_search == TRUE) {
        failures <- c(failures, x_backup[i])
        set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
      }
    }
  }
  
  # handling failures ----
  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0 & initial_search == TRUE) {
    options(mo_failures = sort(unique(failures)))
    plural <- c("value", "it", "was")
    if (n_distinct(failures) > 1) {
      plural <- c("values", "them", "were")
    }
    total_failures <- length(x_input[as.character(x_input) %in% as.character(failures) & !x_input %in% c(NA, NULL, NaN)])
    total_n <- length(x_input[!x_input %in% c(NA, NULL, NaN)])
    msg <- paste0(nr2char(n_distinct(failures)), " unique ", plural[1],
                  " (covering ", percent(total_failures / total_n, round = 1, force_zero = TRUE),
                  ") could not be coerced and ", plural[3], " considered 'unknown'")
    if (n_distinct(failures) <= 10) {
      msg <- paste0(msg, ": ", paste('"', unique(failures), '"', sep = "", collapse = ', '))
    }
    msg <- paste0(msg,  ". Use mo_failures() to review ", plural[2], ". Edit the `allow_uncertain` parameter if needed (see ?as.mo).")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  # handling uncertainties ----
  if (NROW(uncertainties) > 0 & initial_search == TRUE) {
    options(mo_uncertainties = as.list(distinct(uncertainties, input, .keep_all = TRUE)))
    
    plural <- c("", "it")
    if (NROW(uncertainties) > 1) {
      plural <- c("s", "them")
    }
    msg <- paste0("\nResult", plural[1], " of ", nr2char(NROW(uncertainties)), " value", plural[1],
                  " was guessed with uncertainty. Use mo_uncertainties() to review ", plural[2], ".")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  
  # Becker ----
  if (Becker == TRUE | Becker == "all") {
    # See Source. It's this figure:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4187637/figure/F3/
    MOs_staph <- microorganismsDT[genus == "Staphylococcus"]
    setkey(MOs_staph, species)
    CoNS <- MOs_staph[species %in% c("arlettae", "auricularis", "capitis",
                                     "caprae", "carnosus", "chromogenes", "cohnii", "condimenti",
                                     "devriesei", "epidermidis", "equorum", "felis",
                                     "fleurettii", "gallinarum", "haemolyticus",
                                     "hominis", "jettensis", "kloosii", "lentus",
                                     "lugdunensis", "massiliensis", "microti",
                                     "muscae", "nepalensis", "pasteuri", "petrasii",
                                     "pettenkoferi", "piscifermentans", "rostri",
                                     "saccharolyticus", "saprophyticus", "sciuri",
                                     "stepanovicii", "simulans", "succinus",
                                     "vitulinus", "warneri", "xylosus")
                      | (species == "schleiferi" & subspecies %in% c("schleiferi", "")), ..property][[1]]
    CoPS <- MOs_staph[species %in% c("simiae", "agnetis",
                                     "delphini", "lutrae",
                                     "hyicus", "intermedius",
                                     "pseudintermedius", "pseudointermedius",
                                     "schweitzeri", "argenteus")
                      | (species == "schleiferi" & subspecies == "coagulans"), ..property][[1]]
    
    # warn when species found that are not in Becker (2014, PMID 25278577) and Becker (2019, PMID 30872103)
    post_Becker <- c("argensis", "caeli", "cornubiensis", "edaphicus")
    if (any(x %in% MOs_staph[species %in% post_Becker, ..property][[1]])) {
      
      warning("Becker ", italic("et al."), " (2014, 2019) does not contain these species named after their publication: ",
              italic(paste("S.",
                           sort(mo_species(unique(x[x %in% MOs_staph[species %in% post_Becker, ..property][[1]]]))),
                           collapse = ", ")),
              ".",
              call. = FALSE,
              immediate. = TRUE)
    }
    
    x[x %in% CoNS] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
    x[x %in% CoPS] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
    if (Becker == "all") {
      x[x %in% microorganismsDT[mo %like% '^B_STPHY_AUR', ..property][[1]]] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
    }
  }
  
  # Lancefield ----
  if (Lancefield == TRUE | Lancefield == "all") {
    # group A - S. pyogenes
    x[x == microorganismsDT[mo == 'B_STRPT_PYO', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRA', ..property][[1]][1L]
    # group B - S. agalactiae
    x[x == microorganismsDT[mo == 'B_STRPT_AGA', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRB', ..property][[1]][1L]
    # group C
    S_groupC <- microorganismsDT %>% filter(genus == "Streptococcus",
                                            species %in% c("equisimilis", "equi",
                                                           "zooepidemicus", "dysgalactiae")) %>%
      pull(property)
    x[x %in% S_groupC] <- microorganismsDT[mo == 'B_STRPT_GRC', ..property][[1]][1L]
    if (Lancefield == "all") {
      # all Enterococci
      x[x %like% "^(Enterococcus|B_ENTRC)"] <- microorganismsDT[mo == 'B_STRPT_GRD', ..property][[1]][1L]
    }
    # group F - S. anginosus
    x[x == microorganismsDT[mo == 'B_STRPT_ANG', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRF', ..property][[1]][1L]
    # group H - S. sanguinis
    x[x == microorganismsDT[mo == 'B_STRPT_SAN', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRH', ..property][[1]][1L]
    # group K - S. salivarius
    x[x == microorganismsDT[mo == 'B_STRPT_SAL', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRK', ..property][[1]][1L]
  }
  
  # Wrap up ----------------------------------------------------------------
  
  # comply to x, which is also unique and without empty values
  x_input_unique_nonempty <- unique(x_input[!is.na(x_input)
                                            & !is.null(x_input)
                                            & !identical(x_input, "")
                                            & !identical(x_input, "xxx")])
  
  # left join the found results to the original input values (x_input)
  df_found <- data.frame(input = as.character(x_input_unique_nonempty),
                         found = as.character(x),
                         stringsAsFactors = FALSE)
  df_input <- data.frame(input = as.character(x_input),
                         stringsAsFactors = FALSE)
  
  suppressWarnings(
    x <- df_input %>%
      left_join(df_found,
                by = "input") %>%
      pull(found)
  )
  
  if (property == "mo") {
    x <- to_class_mo(x)
  }
  
  if (length(mo_renamed()) > 0) {
    print(mo_renamed())
  }
  
  x
}

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

#' Transform to microorganism ID
#'
#' Use this function to determine a valid microorganism ID ([`mo`]). Determination is done using intelligent rules and the complete taxonomic kingdoms Bacteria, Chromista, Protozoa, Archaea and most microbial species from the kingdom Fungi (see Source). The input can be almost anything: a full name (like `"Staphylococcus aureus"`), an abbreviated name (like `"S. aureus"`), an abbreviation known in the field (like `"MRSA"`), or just a genus. Please see *Examples*.
#' @inheritSection lifecycle Stable lifecycle
#' @param x a character vector or a [`data.frame`] with one or two columns
#' @param Becker a logical to indicate whether *Staphylococci* should be categorised into coagulase-negative *Staphylococci* ("CoNS") and coagulase-positive *Staphylococci* ("CoPS") instead of their own species, according to Karsten Becker *et al.* (1,2). Note that this does not include species that were newly named after these publications, like *S. caeli*.
#'
#' This excludes *Staphylococcus aureus* at default, use `Becker = "all"` to also categorise *S. aureus* as "CoPS".
#' @param Lancefield a logical to indicate whether beta-haemolytic *Streptococci* should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield (3). These *Streptococci* will be categorised in their first group, e.g. *Streptococcus dysgalactiae* will be group C, although officially it was also categorised into groups G and L.
#'
#' This excludes *Enterococci* at default (who are in group D), use `Lancefield = "all"` to also categorise all *Enterococci* as group D.
#' @param allow_uncertain a number between `0` (or `"none"`) and `3` (or `"all"`), or `TRUE` (= `2`) or `FALSE` (= `0`) to indicate whether the input should be checked for less probable results, please see *Details*
#' @param reference_df a [`data.frame`] to be used for extra reference when translating `x` to a valid [`mo`]. See [set_mo_source()] and [get_mo_source()] to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ... other parameters passed on to functions
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
#' ## General info
#' 
#' A microorganism ID from this package (class: [`mo`]) typically looks like these examples:
#' ```
#'   Code               Full name
#'   ---------------    --------------------------------------
#'   B_KLBSL            Klebsiella
#'   B_KLBSL_PNMN       Klebsiella pneumoniae
#'   B_KLBSL_PNMN_RHNS  Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |    |
#'   |   |    |    |
#'   |   |    |     ---> subspecies, a 4-5 letter acronym
#'   |   |     ----> species, a 4-5 letter acronym
#'   |    ----> genus, a 5-7 letter acronym
#'    ----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                             C (Chromista), F (Fungi), P (Protozoa)
#' ```
#'
#' Values that cannot be coered will be considered 'unknown' and will get the MO code `UNKNOWN`.
#'
#' Use the [`mo_*`][mo_property()] functions to get properties based on the returned code, see Examples.
#'
#' The algorithm uses data from the Catalogue of Life (see below) and from one other source (see [microorganisms]).
#'
#' The [as.mo()] function uses several coercion rules for fast and logical results. It assesses the input matching criteria in the following order:
#' 
#' 1. Human pathogenic prevalence: the function  starts with more prevalent microorganisms, followed by less prevalent ones;
#' 2. Taxonomic kingdom: the function starts with determining Bacteria, then Fungi, then Protozoa, then others;
#' 3. Breakdown of input values to identify possible matches.
#'
#' This will lead to the effect that e.g. `"E. coli"` (a microorganism highly prevalent in humans) will return the microbial ID of *Escherichia coli* and not *Entamoeba coli* (a microorganism less prevalent in humans), although the latter would alphabetically come first. 
#' 
#' ## Coping with uncertain results
#' 
#' In addition, the [as.mo()] function can differentiate four levels of uncertainty to guess valid results: 
#' - Uncertainty level 0: no additional rules are applied;
#' - Uncertainty level 1: allow previously accepted (but now invalid) taxonomic names and minor spelling errors;
#' - Uncertainty level 2: allow all of level 1, strip values between brackets, inverse the words of the input, strip off text elements from the end keeping at least two elements;
#' - Uncertainty level 3: allow all of level 1 and 2, strip off text elements from the end, allow any part of a taxonomic name.
#' 
#' This leads to e.g.:
#' - `"Streptococcus group B (known as S. agalactiae)"`. The text between brackets will be removed and a warning will be thrown that the result *Streptococcus group B* (``r as.mo("Streptococcus group B")``) needs review.
#' - `"S. aureus - please mind: MRSA"`. The last word will be stripped, after which the function will try to find a match. If it does not, the second last word will be stripped, etc. Again, a warning will be thrown that the result *Staphylococcus aureus* (``r as.mo("Staphylococcus aureus")``) needs review.
#' - `"Fluoroquinolone-resistant Neisseria gonorrhoeae"`. The first word will be stripped, after which the function will try to find a match. A warning will be thrown that the result *Neisseria gonorrhoeae* (``r as.mo("Neisseria gonorrhoeae")``) needs review.
#'
#' The level of uncertainty can be set using the argument `allow_uncertain`. The default is `allow_uncertain = TRUE`, which is equal to uncertainty level 2. Using `allow_uncertain = FALSE` is equal to uncertainty level 0 and will skip all rules. You can also use e.g. `as.mo(..., allow_uncertain = 1)` to only allow up to level 1 uncertainty.
#' 
#' There are three helper functions that can be run after then [as.mo()] function:
#' - Use [mo_uncertainties()] to get a [`data.frame`] with all values that were coerced to a valid value, but with uncertainty. The output contains a score, that is calculated as \eqn{(n - 0.5 * L) / n}, where *n* is the number of characters of the returned full name of the microorganism, and *L* is the [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) between that full name and the user input.
#' - Use [mo_failures()] to get a [`vector`] with all values that could not be coerced to a valid value.
#' - Use [mo_renamed()] to get a [`data.frame`] with all values that could be coerced based on an old, previously accepted taxonomic name.
#'
#' ## Microbial prevalence of pathogens in humans
#' 
#' The intelligent rules consider the prevalence of microorganisms in humans grouped into three groups, which is available as the `prevalence` columns in the [microorganisms] and [microorganisms.old] data sets. The grouping into prevalence groups is based on experience from several microbiological laboratories in the Netherlands in conjunction with international reports on pathogen prevalence.
#' 
#' Group 1 (most prevalent microorganisms) consists of all microorganisms where the taxonomic class is Gammaproteobacteria or where the taxonomic genus is  *Enterococcus*, *Staphylococcus* or *Streptococcus*. This group consequently contains all common Gram-negative bacteria, such as *Pseudomonas* and *Legionella* and all species within the order Enterobacteriales. 
#' 
#' Group 2 consists of all microorganisms where the taxonomic phylum is Proteobacteria, Firmicutes, Actinobacteria or Sarcomastigophora, or where the taxonomic genus is *Aspergillus*, *Bacteroides*, *Candida*, *Capnocytophaga*, *Chryseobacterium*, *Cryptococcus*, *Elisabethkingia*, *Flavobacterium*, *Fusobacterium*, *Giardia*, *Leptotrichia*, *Mycoplasma*, *Prevotella*, *Rhodotorula*, *Treponema*, *Trichophyton* or *Ureaplasma*. 
#' 
#' Group 3 (least prevalent microorganisms) consists of all other microorganisms.
#' @inheritSection catalogue_of_life Catalogue of Life
#  (source as a section here, so it can be inherited by other man pages:)
#' @section Source:
#' 1. Becker K *et al.* **Coagulase-Negative Staphylococci**. 2014. Clin Microbiol Rev. 27(4): 870–926. <https://dx.doi.org/10.1128/CMR.00109-13>
#' 2. Becker K *et al.* **Implications of identifying the recently defined members of the *S. aureus* complex, *S. argenteus* and *S. schweitzeri*: A position paper of members of the ESCMID Study Group for staphylococci and Staphylococcal Diseases (ESGS).** 2019. Clin Microbiol Infect. <https://doi.org/10.1016/j.cmi.2019.02.028>
#' 3. Lancefield RC **A serological differentiation of human and other groups of hemolytic streptococci**. 1933. J Exp Med. 57(4): 571–95. <https://dx.doi.org/10.1084/jem.57.4.571>
#' 4. Catalogue of Life: Annual Checklist (public online taxonomic database), <http://www.catalogueoflife.org> (check included annual version with [catalogue_of_life_version()]).
#' @export
#' @return A [`character`] vector with class [`mo`]
#' @seealso [microorganisms] for the [`data.frame`] that is being used to determine ID's.
#' 
#' The [mo_property()] functions (like [mo_genus()], [mo_gramstain()]) to get properties based on the returned code.
#' @inheritSection AMR Read more on our website!
#' @examples
#' \donttest{
#' # These examples all return "B_STPHY_AURS", the ID of S. aureus:
#' as.mo("sau") # WHONET code
#' as.mo("stau")
#' as.mo("STAU")
#' as.mo("staaur")
#' as.mo("S. aureus")
#' as.mo("S aureus")
#' as.mo("Staphylococcus aureus")
#' as.mo("Staphylococcus aureus (MRSA)")
#' as.mo("Zthafilokkoockus oureuz") # handles incorrect spelling
#' as.mo("MRSA")    # Methicillin Resistant S. aureus
#' as.mo("VISA")    # Vancomycin Intermediate S. aureus
#' as.mo("VRSA")    # Vancomycin Resistant S. aureus
#' as.mo(115329001) # SNOMED CT code
#' 
#' # Dyslexia is no problem - these all work:
#' as.mo("Ureaplasma urealyticum")
#' as.mo("Ureaplasma urealyticus")
#' as.mo("Ureaplasmium urealytica")
#' as.mo("Ureaplazma urealitycium")
#'
#' as.mo("Streptococcus group A")
#' as.mo("GAS") # Group A Streptococci
#' as.mo("GBS") # Group B Streptococci
#'
#' as.mo("S. epidermidis")                 # will remain species: B_STPHY_EPDR
#' as.mo("S. epidermidis", Becker = TRUE)  # will not remain species: B_STPHY_CONS
#'
#' as.mo("S. pyogenes")                    # will remain species: B_STRPT_PYGN
#' as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRPA
#'
#' # All mo_* functions use as.mo() internally too (see ?mo_property):
#' mo_genus("E. coli")           # returns "Escherichia"
#' mo_gramstain("E. coli")       # returns "Gram negative"
#' 
#' }
#' \dontrun{
#' df$mo <- as.mo(df$microorganism_name)
#'
#' # the select function of tidyverse is also supported:
#' library(dplyr)
#' df$mo <- df %>%
#'   select(microorganism_name) %>%
#'   as.mo()
#'
#' # and can even contain 2 columns, which is convenient for genus/species combinations:
#' df$mo <- df %>%
#'   select(genus, species) %>%
#'   as.mo()
#' # although this works easier and does the same:
#' df <- df %>%
#'   mutate(mo = as.mo(paste(genus, species)))
#' }
as.mo <- function(x, 
                  Becker = FALSE, 
                  Lancefield = FALSE, 
                  allow_uncertain = TRUE, 
                  reference_df = get_mo_source(), 
                  ...) {
  
  check_dataset_integrity()
  
  # start off with replaced language-specific non-ASCII characters with ASCII characters
  x <- parse_and_convert(x)
  
  # WHONET: xxx = no growth
  x[tolower(as.character(paste0(x, ""))) %in% c("", "xxx", "na", "nan")] <- NA_character_
  # Laboratory systems: remove entries like "no growth" etc
  x[trimws2(x) %like% "(no .*growth|keine? .*wachtstum|geen .*groei|no .*crecimientonon|sem .*crescimento|pas .*croissance)"] <- NA_character_
  x[trimws2(x) %like% "^(no|not|kein|geen|niet|non|sem) [a-z]+"] <- "UNKNOWN"
  
  uncertainty_level <- translate_allow_uncertain(allow_uncertain)
  
  if (mo_source_isvalid(reference_df)
      & isFALSE(Becker)
      & isFALSE(Lancefield)
      & !is.null(reference_df)
      & all(x %in% reference_df[, 1][[1]])) {
    
    # has valid own reference_df
    # (data.table not faster here)
    reference_df <- reference_df %>% filter(!is.na(mo))
    # keep only first two columns, second must be mo
    if (colnames(reference_df)[1] == "mo") {
      reference_df <- reference_df[, c(2, 1)]
    } else {
      reference_df <- reference_df[, c(1, 2)]
    }
    colnames(reference_df)[1] <- "x"
    # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
    suppressWarnings(
      y <- data.frame(x = x, stringsAsFactors = FALSE) %>%
        left_join(reference_df, by = "x") %>%
        pull("mo")
    )
    
  } else if (all(x %in% MO_lookup$mo)
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    y <- x
    
  } else {
    # will be checked for mo class in validation and uses exec_as.mo internally if necessary
    y <- mo_validate(x = x, property = "mo",
                     Becker = Becker, Lancefield = Lancefield,
                     allow_uncertain = uncertainty_level, reference_df = reference_df,
                     ...)
  }
  
  to_class_mo(y)
}

to_class_mo <- function(x) {
  structure(.Data = x,
            class = c("mo", "character"))
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  inherits(x, "mo")
}

# param property a column name of microorganisms
# param initial_search logical - is FALSE when coming from uncertain tries, which uses exec_as.mo internally too
# param dyslexia_mode logical - also check for characters that resemble others
# param debug logical - show different lookup texts while searching
# param reference_data_to_use data.frame - the data set to check for
exec_as.mo <- function(x,
                       Becker = FALSE,
                       Lancefield = FALSE,
                       allow_uncertain = TRUE,
                       reference_df = get_mo_source(),
                       property = "mo",
                       initial_search = TRUE,
                       dyslexia_mode = FALSE,
                       debug = FALSE,
                       reference_data_to_use = MO_lookup) {
  
  check_dataset_integrity()
  
  lookup <- function(needle, column = property, haystack = reference_data_to_use, n = 1, debug_mode = debug) {
    # `column` can be NULL for all columns, or a selection
    # returns a character (vector) - if `column` > length 1 then with columns as names
    if (isTRUE(debug_mode)) {
      cat(font_silver("looking up: ", substitute(needle), "\n", collapse = ""))
    }
    if (length(column) == 1) {
      res <- haystack[which(eval(substitute(needle), envir = haystack, enclos = parent.frame())), column, drop = TRUE]
      res <- as.character(res)
      if (length(res) == 0) {
        NA_character_
      } else {
        res[seq_len(min(n, length(res)))]
      }
    } else {
      if (is.null(column)) {
        column <- names(haystack)
      }
      res <- haystack[which(eval(substitute(needle), envir = haystack, enclos = parent.frame())), , drop = FALSE]
      res <- res[seq_len(min(n, nrow(res))), column, drop = TRUE]
      if (NROW(res) == 0) {
        res <- rep(NA_character_, length(column))
      }
      res <- as.character(res)
      names(res) <- column
      res
    }
  }
  
  # start off with replaced language-specific non-ASCII characters with ASCII characters
  x <- parse_and_convert(x)
  
  # WHONET: xxx = no growth
  x[tolower(as.character(paste0(x, ""))) %in% c("", "xxx", "na", "nan")] <- NA_character_
  # Laboratory systems: remove entries like "no growth" etc
  x[trimws2(x) %like% "(no .*growth|keine? .*wachtstum|geen .*groei|no .*crecimientonon|sem .*crescimento|pas .*croissance)"] <- NA_character_
  x[trimws2(x) %like% "^(no|not|kein|geen|niet|non|sem) [a-z]+"] <- "UNKNOWN"
  
  if (initial_search == TRUE) {
    options(mo_failures = NULL)
    options(mo_uncertainties = NULL)
    options(mo_renamed = NULL)
  }
  options(mo_renamed_last_run = NULL)
  
  uncertainties <- data.frame(uncertainty = integer(0),
                              input = character(0),
                              fullname = character(0),
                              renamed_to = character(0),
                              mo = character(0), 
                              stringsAsFactors = FALSE)
  failures <- character(0)
  uncertainty_level <- translate_allow_uncertain(allow_uncertain)
  old_mo_warning <- FALSE
  
  x_input <- x
  # already strip leading and trailing spaces
  x <- trimws(x)
  # only check the uniques, which is way faster
  x <- unique(x)
  # remove empty values (to later fill them in again with NAs)
  # ("xxx" is WHONET code for 'no growth')
  x <- x[!is.na(x)
         & !is.null(x)
         & !identical(x, "")
         & !identical(x, "xxx")]
  
  # conversion of old MO codes from v0.5.0 (ITIS) to later versions (Catalogue of Life)
  if (any(x %like_case% "^[BFP]_[A-Z]{3,7}") & !all(x %in% microorganisms$mo)) {
    leftpart <- gsub("^([BFP]_[A-Z]{3,7}).*", "\\1", x)
    if (any(leftpart %in% names(mo_codes_v0.5.0))) {
      old_mo_warning <- TRUE
      rightpart <- gsub("^[BFP]_[A-Z]{3,7}(.*)", "\\1", x)
      leftpart <- mo_codes_v0.5.0[leftpart]
      x[!is.na(leftpart)] <- paste0(leftpart[!is.na(leftpart)], rightpart[!is.na(leftpart)])
    }
    # now check if some are still old
    still_old <- x[x %in% names(mo_codes_v0.5.0)]
    if (length(still_old) > 0) {
      old_mo_warning <- TRUE
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
  
  # defined df to check for
  if (!is.null(reference_df)) {
    if (!mo_source_isvalid(reference_df)) {
      stop("`reference_df` must contain a column `mo` with values from the 'microorganisms' data set.", call. = FALSE)
    }
    reference_df <- reference_df %>% filter(!is.na(mo))
    # keep only first two columns, second must be mo
    if (colnames(reference_df)[1] == "mo") {
      reference_df <- reference_df[, c(2, 1)]
    } else {
      reference_df <- reference_df[, c(1, 2)]
    }
    colnames(reference_df)[1] <- "x"
    # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
  }
  
  # all empty
  if (all(identical(trimws(x_input), "") | is.na(x_input) | length(x) == 0)) {
    if (property == "mo") {
      return(to_class_mo(rep(NA_character_, length(x_input))))
    } else {
      return(rep(NA_character_, length(x_input)))
    }
    
  } else if (all(x %in% reference_df[, 1][[1]])) {
    # all in reference df
    colnames(reference_df)[1] <- "x"
    suppressWarnings(
      x <- data.frame(x = x, stringsAsFactors = FALSE) %>%
        left_join(reference_df, by = "x") %>%
        left_join(microorganisms, by = "mo") %>%
        pull(property)
    )
    
  } else if (all(x %in% reference_data_to_use$mo)) {
    x <- data.frame(mo = x, stringsAsFactors = FALSE) %>% 
      left_join_microorganisms(by = "mo") %>% 
      pull(property)
    
  } else if (all(tolower(x) %in% reference_data_to_use$fullname_lower)) {
    # we need special treatment for very prevalent full names, they are likely!
    # e.g. as.mo("Staphylococcus aureus")
    x <- data.frame(fullname_lower = tolower(x), stringsAsFactors = FALSE) %>% 
      left_join_MO_lookup(by = "fullname_lower") %>% 
      pull(property)
    # x <- reference_data_to_use[data.table(fullname_lower = tolower(x)),
    #                            on = "fullname_lower",
    #                            ..property][[1]]
    
  } else if (all(toupper(x) %in% microorganisms.codes$code)) {
    # commonly used MO codes
    x <- data.frame(code = toupper(x), stringsAsFactors = FALSE) %>%
      left_join(microorganisms.codes, by = "code") %>%
      left_join_MO_lookup(by = "mo") %>%
      pull(property)
    # y <- as.data.table(microorganisms.codes)[data.table(code = toupper(x)),
    #                                               on = "code", ]
    # 
    # x <- reference_data_to_use[data.table(mo = y[["mo"]]),
    #                            on = "mo",
    #                            ..property][[1]]
    
  } else if (all(x %in% microorganisms.translation$mo_old)) {
    # is an old mo code, used in previous versions of this package
    old_mo_warning <- TRUE
    x <- data.frame(mo_old = toupper(x), stringsAsFactors = FALSE) %>%
      left_join(microorganisms.translation, by = "mo_old") %>%
      rename(mo = mo_new) %>% 
      left_join_MO_lookup(by = "mo") %>%
      pull(property)
    
  } else if (!all(x %in% microorganisms[, property])) {
    
    strip_whitespace <- function(x, dyslexia_mode) {
      # all whitespaces (tab, new lines, etc.) should be one space
      # and spaces before and after should be omitted
      trimmed <- trimws2(x)
      # also, make sure the trailing and leading characters are a-z or 0-9
      # in case of non-regex
      if (dyslexia_mode == FALSE) {
        trimmed <- gsub("^[^a-zA-Z0-9)(]+", "", trimmed)
        trimmed <- gsub("[^a-zA-Z0-9)(]+$", "", trimmed)
      }
      trimmed
    }
    
    x_backup_untouched <- x
    x <- strip_whitespace(x, dyslexia_mode)
    x_backup <- x
    
    # from here on case-insensitive
    x <- tolower(x)
    
    x_backup[grepl("^(fungus|fungi)$", x)] <- "F_FUNGUS" # will otherwise become the kingdom
    
    # remove spp and species
    x <- gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x)
    x <- gsub("(spp.?|subsp.?|subspecies|biovar|serovar|species)", "", x)
    x <- gsub("^([a-z]{2,4})(spe.?)$", "\\1", x) # when ending in SPE instead of SPP and preceded by 2-4 characters
    x <- strip_whitespace(x, dyslexia_mode)
    
    x_backup_without_spp <- x
    x_species <- paste(x, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(gruppe|groep|grupo|gruppo|groupe)", "group", x)
    # no groups and complexes as ending
    x <- gsub("(complex|group)$", "", x)
    x <- gsub("((an)?aero+b)[a-z]*", "", x)
    x <- gsub("^atyp[a-z]*", "", x)
    x <- gsub("(vergroen)[a-z]*", "viridans", x)
    x <- gsub("[a-z]*diff?erent[a-z]*", "", x)
    x <- gsub("(hefe|gist|gisten|levadura|lievito|fermento|levure)[a-z]*", "yeast", x)
    x <- gsub("(schimmels?|mofo|molde|stampo|moisissure|fungi)[a-z]*", "fungus", x)
    x <- gsub("fungus[ph|f]rya", "fungiphrya", x)
    # no contamination
    x <- gsub("(contamination|kontamination|mengflora|contaminaci.n|contamina..o)", "", x)
    # remove non-text in case of "E. coli" except dots and spaces
    x <- trimws(gsub("[^.a-zA-Z0-9/ \\-]+", " ", x))
    # but make sure that dots are followed by a space
    x <- gsub("[.] ?", ". ", x)
    # replace minus by a space
    x <- gsub("-+", " ", x)
    # replace hemolytic by haemolytic
    x <- gsub("ha?emoly", "haemoly", x)
    # place minus back in streptococci
    x <- gsub("(alpha|beta|gamma).?ha?emoly", "\\1-haemoly", x)
    # remove genus as first word
    x <- gsub("^genus ", "", x)
    # remove 'uncertain'-like texts
    x <- trimws(gsub("(uncertain|susp[ie]c[a-z]+|verdacht)", "", x))
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
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+a\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+u\\+m)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("e+", "e+", x)
      x <- gsub("o+", "o+", x)
      x <- gsub("(.)\\1+", "\\1+", x)
      # allow multiplication of all other consonants
      x <- gsub("([bdgjlnrw]+)", "\\1+", x)
      # allow ending in -en or -us
      x <- gsub("e\\+n(?![a-z[])", "(e+n|u+(c|k|q|qu|s|z|x|ks)+)", x, perl = TRUE)
      # if the input is longer than 10 characters, allow any forgotten consonant between all characters, as some might just have forgotten one...
      # this will allow "Pasteurella damatis" to be correctly read as "Pasteurella dagmatis".
      consonants <- paste(letters[!letters %in% c("a", "e", "i", "o", "u")], collapse = "")
      x[nchar(x_backup_without_spp) > 10] <- gsub("[+]", paste0("+[", consonants, "]?"), x[nchar(x_backup_without_spp) > 10])
      # allow au and ou after all these regex implementations
      x <- gsub("a+[bcdfghjklmnpqrstvwxyz]?u+[bcdfghjklmnpqrstvwxyz]?", "(a+u+|o+u+)[bcdfghjklmnpqrstvwxyz]?", x, fixed = TRUE)
      x <- gsub("o+[bcdfghjklmnpqrstvwxyz]?u+[bcdfghjklmnpqrstvwxyz]?", "(a+u+|o+u+)[bcdfghjklmnpqrstvwxyz]?", x, fixed = TRUE)
    }
    x <- strip_whitespace(x, dyslexia_mode)
    # make sure to remove regex overkill (will lead to errors)
    x <- gsub("++", "+", x, fixed = TRUE)
    x <- gsub("?+", "?", x, fixed = TRUE)
    
    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    x_trimmed_without_group <- gsub(" gro.u.p$", "", x_trimmed)
    # remove last part from "-" or "/"
    x_trimmed_without_group <- gsub("(.*)[-/].*", "\\1", x_trimmed_without_group)
    # replace space and dot by regex sign
    x_withspaces <- gsub("[ .]+", ".* ", x)
    x <- gsub("[ .]+", ".*", x)
    # add start en stop regex
    x <- paste0("^", x, "$")
    
    x_withspaces_start_only <- paste0("^", x_withspaces)
    x_withspaces_end_only <- paste0(x_withspaces, "$")
    x_withspaces_start_end <- paste0("^", x_withspaces, "$")
    
    if (isTRUE(debug)) {
      cat(paste0(font_blue("x"), '                       "', x, '"\n'))
      cat(paste0(font_blue("x_species"), '               "', x_species, '"\n'))
      cat(paste0(font_blue("x_withspaces_start_only"), ' "', x_withspaces_start_only, '"\n'))
      cat(paste0(font_blue("x_withspaces_end_only"), '   "', x_withspaces_end_only, '"\n'))
      cat(paste0(font_blue("x_withspaces_start_end"), '  "', x_withspaces_start_end, '"\n'))
      cat(paste0(font_blue("x_backup"), '                "', x_backup, '"\n'))
      cat(paste0(font_blue("x_backup_without_spp"), '    "', x_backup_without_spp, '"\n'))
      cat(paste0(font_blue("x_trimmed"), '               "', x_trimmed, '"\n'))
      cat(paste0(font_blue("x_trimmed_species"), '       "', x_trimmed_species, '"\n'))
      cat(paste0(font_blue("x_trimmed_without_group"), ' "', x_trimmed_without_group, '"\n'))
    }
    
    if (initial_search == TRUE) {
      progress <- progress_estimated(n = length(x), n_min = 25) # start if n >= 25
      on.exit(close(progress))
    }
    
    for (i in seq_len(length(x))) {
      
      if (initial_search == TRUE) {
        progress$tick()
      }

      # valid MO code ----
      found <- lookup(mo == toupper(x_backup[i]))
      if (!is.na(found)) {
        x[i] <- found[1L]
        next
      }
      
      # valid fullname ----
      found <- lookup(fullname_lower %in% gsub("[^a-zA-Z0-9_. -]", "", tolower(c(x_backup[i], x_backup_without_spp[i]))))
      # added the gsub() for "(unknown fungus)", since fullname_lower does not contain brackets
      if (!is.na(found)) {
        x[i] <- found[1L]
        next
      }
      
      # old fullname ----
      found <- lookup(fullname_lower %in% tolower(c(x_backup[i], x_backup_without_spp[i])),
                      column = NULL, # all columns
                      haystack = MO.old_lookup)
      if (!all(is.na(found))) {
        # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
        # mo_ref() of "Chlamydia psittaci" will be "Page, 1968" (with warning)
        # mo_ref() of "Chlamydophila psittaci" will be "Everett et al., 1999"
        if (property == "ref") {
          x[i] <- found["ref"]
        } else {
          x[i] <- lookup(fullname == found["fullname_new"], haystack = MO_lookup)
        }
        options(mo_renamed_last_run = found["fullname"])
        was_renamed(name_old = found["fullname"],
                    name_new = lookup(fullname == found["fullname_new"], "fullname", haystack = MO_lookup),
                    ref_old = found["ref"],
                    ref_new = lookup(fullname == found["fullname_new"], "ref", haystack = MO_lookup),
                    mo = lookup(fullname == found["fullname_new"], "mo", haystack = MO_lookup))
        next
      }
      
      # old mo code, used in previous versions of this package ----
      if (x_backup[i] %in% microorganisms.translation$mo_old) {
        old_mo_warning <- TRUE
        found <- lookup(mo_old == toupper(x_backup[i]), column = "mo_new", haystack = microorganisms.translation)
        found <- lookup(mo == found)
        if (!is.na(found)) {
          # get property
          x[i] <- found[1L]
          next
        }
      }
      
      if (x_backup[i] %like_case% "\\(unknown [a-z]+\\)" | tolower(x_backup_without_spp[i]) %in% c("other", "none", "unknown")) {
        # empty and nonsense values, ignore without warning
        x[i] <- lookup(mo == "UNKNOWN")
        next
      }
      
      # exact SNOMED code ----
      if (x_backup[i] %like% "^[0-9]+$") {
        snomed_found <- unlist(lapply(reference_data_to_use$snomed,
                                      function(s) if (x_backup[i] %in% s) {
                                        TRUE
                                      } else {
                                        FALSE
                                      }))
        if (sum(snomed_found, na.rm = TRUE) > 0) {
          found <- reference_data_to_use[snomed_found == TRUE, property][[1]]
          if (!is.na(found)) {
            x[i] <- found[1L]
            next
          }
        }
      }
      
      # very probable: is G. species ----
      found <- lookup(g_species %in% gsub("[^a-z0-9/ \\-]+", "", 
                                          tolower(c(x_backup[i], x_backup_without_spp[i]))))
      if (!is.na(found)) {
        x[i] <- found[1L]
        next
      }
      
      # WHONET and other common LIS codes ----
      found <- lookup(code %in% toupper(c(x_backup_untouched[i], x_backup[i], x_backup_without_spp[i])),
                      column = "mo", 
                      haystack = microorganisms.codes)
      if (!is.na(found)) {
        x[i] <- lookup(mo == found)
        next
      }
      
      # user-defined reference ----
      if (!is.null(reference_df)) {
        if (x_backup[i] %in% reference_df[, 1]) {
          # already checked integrity of reference_df, all MOs are valid
          ref_mo <- reference_df[reference_df[, 1] == x_backup[i], "mo"][[1L]]
          x[i] <- lookup(mo == ref_mo)
          next
        }
      }
      
      # WHONET: xxx = no growth
      if (tolower(as.character(paste0(x_backup_without_spp[i], ""))) %in% c("", "xxx", "na", "nan")) {
        x[i] <- NA_character_
        next
      }
      
      # check for very small input, but ignore the O antigens of E. coli
      if (nchar(gsub("[^a-zA-Z]", "", x_trimmed[i])) < 3
          & !toupper(x_backup_without_spp[i]) %like_case% "O?(26|103|104|104|111|121|145|157)") {
        # fewer than 3 chars and not looked for species, add as failure
        x[i] <- lookup(mo == "UNKNOWN")
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
        }
        next
      }
      
      if (x_backup_without_spp[i] %like_case% "(virus|viridae)") {
        # there is no fullname like virus or viridae, so don't try to coerce it
        x[i] <- NA_character_
        next
      }
      
      # translate known trivial abbreviations to genus + species ----
      if (toupper(x_backup_without_spp[i]) %in% c("MRSA", "MSSA", "VISA", "VRSA")
          | x_backup_without_spp[i] %like_case% " (mrsa|mssa|visa|vrsa) ") {
        x[i] <- lookup(fullname == "Staphylococcus aureus")
        next
      }
      if (toupper(x_backup_without_spp[i]) %in% c("MRSE", "MSSE")
          | x_backup_without_spp[i] %like_case% " (mrse|msse) ") {
        x[i] <- lookup(fullname == "Staphylococcus epidermidis")
        next
      }
      if (toupper(x_backup_without_spp[i]) == "VRE"
          | x_backup_without_spp[i] %like_case% " vre "
          | x_backup_without_spp[i] %like_case% "(enterococci|enterokok|enterococo)[a-z]*?$")  {
        x[i] <- lookup(genus == "Enterococcus")
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
          | x_backup_without_spp[i] %like_case% "o?(26|103|104|111|121|145|157)") {
        x[i] <- lookup(fullname == "Escherichia coli")
        next
      }
      if (toupper(x_backup_without_spp[i]) == "MRPA"
          | x_backup_without_spp[i] %like_case% " mrpa ") {
        # multi resistant P. aeruginosa
        x[i] <- lookup(fullname == "Pseudomonas aeruginosa")
        next
      }
      if (toupper(x_backup_without_spp[i]) == "CRSM") {
        # co-trim resistant S. maltophilia
        x[i] <- lookup(fullname == "Stenotrophomonas maltophilia")
        next
      }
      if (toupper(x_backup_without_spp[i]) %in% c("PISP", "PRSP", "VISP", "VRSP")
          | x_backup_without_spp[i] %like_case% " (pisp|prsp|visp|vrsp) ") {
        # peni I, peni R, vanco I, vanco R: S. pneumoniae
        x[i] <- lookup(fullname == "Streptococcus pneumoniae")
        next
      }
      if (x_backup_without_spp[i] %like_case% "^g[abcdfghk]s$") {
        # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRPB)
        x[i] <- lookup(mo == toupper(gsub("g([abcdfghk])s",
                                          "B_STRPT_GRP\\1",
                                          x_backup_without_spp[i])))
        next
      }
      if (x_backup_without_spp[i] %like_case% "(streptococ|streptokok).* [abcdfghk]$") {
        # Streptococci in different languages, like "estreptococos grupo B"
        x[i] <- lookup(mo == toupper(gsub(".*(streptococ|streptokok|estreptococ).* ([abcdfghk])$", 
                                          "B_STRPT_GRP\\2",
                                          x_backup_without_spp[i])))
        next
      }
      if (x_backup_without_spp[i] %like_case% "group [abcdfghk] (streptococ|streptokok|estreptococ)") {
        # Streptococci in different languages, like "Group A Streptococci"
        x[i] <- lookup(mo == toupper(gsub(".*group ([abcdfghk]) (streptococ|streptokok|estreptococ).*", 
                                          "B_STRPT_GRP\\1",
                                          x_backup_without_spp[i])))
        next
      }
      if (x_backup_without_spp[i] %like_case% "haemoly.*strept") {
        # Haemolytic streptococci in different languages
        x[i] <- lookup(mo == "B_STRPT_HAEM")
        next
      }
      # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese) ----
      if (x_backup_without_spp[i] %like_case% "[ck]oagulas[ea] negatie?[vf]"
          | x_trimmed[i] %like_case% "[ck]oagulas[ea] negatie?[vf]"
          | x_backup_without_spp[i] %like_case% "[ck]o?ns[^a-z]?$") {
        # coerce S. coagulase negative
        x[i] <- lookup(mo == "B_STPHY_CONS")
        next
      }
      if (x_backup_without_spp[i] %like_case% "[ck]oagulas[ea] positie?[vf]"
          | x_trimmed[i] %like_case% "[ck]oagulas[ea] positie?[vf]"
          | x_backup_without_spp[i] %like_case% "[ck]o?ps[^a-z]?$") {
        # coerce S. coagulase positive
        x[i] <- lookup(mo == "B_STPHY_COPS")
        next
      }
      # streptococcal groups: milleri and viridans
      if (x_trimmed[i] %like_case% "strepto.* mil+er+i"
          | x_backup_without_spp[i] %like_case% "strepto.* mil+er+i"
          | x_backup_without_spp[i] %like_case% "mgs[^a-z]?$") {
        # Milleri Group Streptococcus (MGS)
        x[i] <- lookup(mo == "B_STRPT_MILL")
        next
      }
      if (x_trimmed[i] %like_case% "strepto.* viridans"
          | x_backup_without_spp[i] %like_case% "strepto.* viridans"
          | x_backup_without_spp[i] %like_case% "vgs[^a-z]?$") {
        # Viridans Group Streptococcus (VGS)
        x[i] <- lookup(mo == "B_STRPT_VIRI")
        next
      }
      if (x_backup_without_spp[i] %like_case% "gram[ -]?neg.*"
          | x_backup_without_spp[i] %like_case% "negatie?[vf]"
          | x_trimmed[i] %like_case% "gram[ -]?neg.*") {
        # coerce Gram negatives
        x[i] <- lookup(mo == "B_GRAMN")
        next
      }
      if (x_backup_without_spp[i] %like_case% "gram[ -]?pos.*"
          | x_backup_without_spp[i] %like_case% "positie?[vf]"
          | x_trimmed[i] %like_case% "gram[ -]?pos.*") {
        # coerce Gram positives
        x[i] <- lookup(mo == "B_GRAMP")
        next
      }
      if (x_backup_without_spp[i] %like_case% "mycoba[ck]teri.[nm]?$") {
        # coerce mycobacteria in multiple languages
        x[i] <- lookup(genus == "Mycobacterium")
        next
      }
      
      if (x_backup_without_spp[i] %like_case% "salmonella [a-z]+ ?.*") {
        if (x_backup_without_spp[i] %like_case% "salmonella group") {
          # Salmonella Group A to Z, just return S. species for now
          x[i] <- lookup(genus == "Salmonella")
          next
        } else if (grepl("[sS]almonella [A-Z][a-z]+ ?.*", x_backup[i], ignore.case = FALSE) &
                   !x_backup[i] %like% "t[iy](ph|f)[iy]") {
          # Salmonella with capital letter species like "Salmonella Goettingen" - they're all S. enterica
          # except for S. typhi, S. paratyphi, S. typhimurium
          x[i] <- lookup(fullname == "Salmonella enterica")
          uncertainties <- rbind(uncertainties,
                                 format_uncertainty_as_df(uncertainty_level = 1,
                                                          input = x_backup[i],
                                                          result_mo = lookup(fullname == "Salmonella enterica", "mo")))
          next
        }
      }
      
      # trivial names known to the field:
      if ("meningococcus" %like_case% x_trimmed[i]) {
        # coerce Neisseria meningitidis
        x[i] <- lookup(fullname == "Neisseria meningitidis")
        next
      }
      if ("gonococcus" %like_case% x_trimmed[i]) {
        # coerce Neisseria gonorrhoeae
        x[i] <- lookup(fullname == "Neisseria gonorrhoeae")
        next
      }
      if ("pneumococcus" %like_case% x_trimmed[i]) {
        # coerce Streptococcus penumoniae
        x[i] <- lookup(fullname == "Streptococcus pneumoniae")
        next
      }
      # }
      
      # NOW RUN THROUGH DIFFERENT PREVALENCE LEVELS
      check_per_prevalence <- function(data_to_check,
                                       data.old_to_check,
                                       a.x_backup,
                                       b.x_trimmed,
                                       c.x_trimmed_without_group,
                                       d.x_withspaces_start_end,
                                       e.x_withspaces_start_only,
                                       f.x_withspaces_end_only,
                                       g.x_backup_without_spp,
                                       h.x_species,
                                       i.x_trimmed_species) {
        
        # FIRST TRY FULLNAMES AND CODES ----
        # if only genus is available, return only genus
        
        if (all(!c(x[i], b.x_trimmed) %like_case% " ")) {
          found <- lookup(fullname_lower %in% c(h.x_species, i.x_trimmed_species),
                          haystack = data_to_check)
          if (!is.na(found)) {
            x[i] <- found[1L]
            return(x[i])
          }
          if (nchar(g.x_backup_without_spp) >= 6) {
            found <- lookup(fullname_lower %like_case% paste0("^", unregex(g.x_backup_without_spp), "[a-z]+"),
                            haystack = data_to_check)
            if (!is.na(found)) {
              x[i] <- found[1L]
              return(x[i])
            }
          }
          # rest of genus only is in allow_uncertain part.
        }
        
        # allow no codes less than 4 characters long, was already checked for WHONET earlier
        if (nchar(g.x_backup_without_spp) < 4) {
          x[i] <- lookup(mo == "UNKNOWN")
          if (initial_search == TRUE) {
            failures <- c(failures, a.x_backup)
          }
          return(x[i])
        }
        
        # try probable: trimmed version of fullname ----
        found <- lookup(fullname_lower %in% tolower(g.x_backup_without_spp),
                        haystack = data_to_check)
        if (!is.na(found)) {
          return(found[1L])
        }
        
        # try any match keeping spaces ----
        found <- lookup(fullname_lower %like_case% d.x_withspaces_start_end,
                        haystack = data_to_check)
        if (!is.na(found) & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not ending with $ ----
        found <- lookup(fullname_lower %like_case% paste0(trimws(e.x_withspaces_start_only), " "),
                        haystack = data_to_check)
        if (!is.na(found)) {
          return(found[1L])
        }
        found <- lookup(fullname_lower %like_case% e.x_withspaces_start_only,
                        haystack = data_to_check)
        if (!is.na(found) & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not start with ^ ----
        found <- lookup(fullname_lower %like_case% paste0(" ", trimws(f.x_withspaces_end_only)),
                        haystack = data_to_check)
        if (!is.na(found)) {
          return(found[1L])
        }
        
        # try a trimmed version
        found <- lookup(fullname_lower %like_case% b.x_trimmed |
                          fullname_lower %like_case% c.x_trimmed_without_group,
                        haystack = data_to_check)
        if (!is.na(found) & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        
        # try splitting of characters in the middle and then find ID ----
        # only when text length is 6 or lower
        # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
        if (nchar(g.x_backup_without_spp) <= 6) {
          x_length <- nchar(g.x_backup_without_spp)
          x_split <- paste0("^",
                            g.x_backup_without_spp %>% substr(1, x_length / 2),
                            ".* ",
                            g.x_backup_without_spp %>% substr((x_length / 2) + 1, x_length))
          found <- lookup(fullname_lower %like_case% x_split,
                          haystack = data_to_check)
          if (!is.na(found)) {
            return(found[1L])
          }
        }
        
        # try fullname without start and without nchar limit of >= 6 ----
        # like "K. pneu rhino" >> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
        found <- lookup(fullname_lower %like_case% e.x_withspaces_start_only,
                        haystack = data_to_check)
        if (!is.na(found)) {
          return(found[1L])
        }
        
        # MISCELLANEOUS ----
        
        # look for old taxonomic names ----
        found <- lookup(fullname_lower %like_case% e.x_withspaces_start_only,
                        column = NULL, # all columns
                        haystack = data.old_to_check)
        if (!all(is.na(found))) {
          # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
          # mo_ref() of "Chlamydia psittaci" will be "Page, 1968" (with warning)
          # mo_ref() of "Chlamydophila psittaci" will be "Everett et al., 1999"
          if (property == "ref") {
            x[i] <- found["ref"]
          } else {
            x[i] <- lookup(fullname == found["fullname_new"], haystack = MO_lookup)
          }
          options(mo_renamed_last_run = found["fullname"])
          was_renamed(name_old = found["fullname"],
                      name_new = lookup(fullname == found["fullname_new"], "fullname", haystack = MO_lookup),
                      ref_old = found["ref"],
                      ref_new = lookup(fullname == found["fullname_new"], "ref", haystack = MO_lookup),
                      mo = lookup(fullname == found["fullname_new"], "mo", haystack = MO_lookup))
          return(x[i])
        }
        
        # check for uncertain results ----
        uncertain_fn <- function(a.x_backup,
                                 b.x_trimmed,
                                 d.x_withspaces_start_end,
                                 e.x_withspaces_start_only,
                                 f.x_withspaces_end_only,
                                 g.x_backup_without_spp,
                                 uncertain.reference_data_to_use) {
          
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
              message("Running '", d.x_withspaces_start_end, "' and '", e.x_withspaces_start_only, "'")
            }
            found <- lookup(fullname_lower %like_case% d.x_withspaces_start_end |
                              fullname_lower %like_case% e.x_withspaces_start_only,
                            column = NULL, # all columns
                            haystack = data.old_to_check)
            if (!all(is.na(found)) & nchar(g.x_backup_without_spp) >= 6) {
              if (property == "ref") {
                # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
                # mo_ref("Chlamydia psittaci) = "Page, 1968" (with warning)
                # mo_ref("Chlamydophila psittaci) = "Everett et al., 1999"
                x <- found["ref"]
              } else {
                x <- lookup(fullname == found["fullname_new"], haystack = MO_lookup)
              }
              was_renamed(name_old = found["fullname"],
                          name_new = lookup(fullname == found["fullname_new"], "fullname", haystack = MO_lookup),
                          ref_old = found["ref"],
                          ref_new = lookup(fullname == found["fullname_new"], "ref", haystack = MO_lookup),
                          mo = lookup(fullname == found["fullname_new"], "mo", haystack = MO_lookup))
              options(mo_renamed_last_run = found["fullname"])
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = lookup(fullname == found["fullname_new"], "mo", haystack = MO_lookup)))
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
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found)) {
              found_result <- found
              found <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
            }
          }
          
          # UNCERTAINTY LEVEL 2 ----
          if (uncertainty_level >= 2) {
            now_checks_for_uncertainty_level <- 2
            
            # (3) look for genus only, part of name ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (3) look for genus only, part of name\n")
            }
            if (nchar(g.x_backup_without_spp) > 4 & !b.x_trimmed %like_case% " ") {
              if (!grepl("^[A-Z][a-z]+", b.x_trimmed, ignore.case = FALSE)) {
                if (isTRUE(debug)) {
                  message("Running '", paste(b.x_trimmed, "species"), "'")
                }
                # not when input is like Genustext, because then Neospora would lead to Actinokineospora
                found <- lookup(fullname_lower %like_case% paste(b.x_trimmed, "species"),
                                haystack = uncertain.reference_data_to_use)
                if (!is.na(found)) {
                  found_result <- found
                  found <- lookup(mo == found)
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result))
                  return(found)
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
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
              found_result <- found
              found <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
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
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
              found_result <- found
              found <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
            }
            
            # (6) try to strip off half an element from end and check the remains ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (6) try to strip off half an element from end and check the remains\n")
            }
            x_strip <- a.x_backup %>% strsplit("[ .]") %>% unlist()
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                lastword <- x_strip[length(x_strip) - i + 1]
                lastword_half <- substr(lastword, 1, as.integer(nchar(lastword) / 2))
                # remove last half of the second term
                x_strip_collapsed <- paste(c(x_strip[seq_len(length(x_strip) - i)], lastword_half), collapse = " ")
                if (nchar(x_strip_collapsed) >= 4 & nchar(lastword_half) > 2) {
                  if (isTRUE(debug)) {
                    message("Running '", x_strip_collapsed, "'")
                  }
                  # first try without dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  if (empty_result(found)) {
                    # then with dyslexia mode
                    found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  }
                  if (!empty_result(found)) {
                    found_result <- found
                    found <- lookup(mo == found)
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result))
                    return(found)
                  }
                }
              }
            }
            # (7) try to strip off one element from end and check the remains ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (7) try to strip off one element from end and check the remains\n")
            }
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                x_strip_collapsed <- paste(x_strip[seq_len(length(x_strip) - i)], collapse = " ")
                if (nchar(x_strip_collapsed) >= 6) {
                  if (isTRUE(debug)) {
                    message("Running '", x_strip_collapsed, "'")
                  }
                  # first try without dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  if (empty_result(found)) {
                    # then with dyslexia mode
                    found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  }
                  if (!empty_result(found)) {
                    found_result <- found
                    found <- lookup(mo == found)
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result))
                    return(found)
                  }
                }
              }
            }
            # (8) check for unknown yeasts/fungi ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (8) check for unknown yeasts/fungi\n")
            }
            if (b.x_trimmed %like_case% "yeast") {
              found <- "F_YEAST"
              found_result <- found
              found <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
            }
            if (b.x_trimmed %like_case% "(fungus|fungi)" & !b.x_trimmed %like_case% "fungiphrya") {
              found <- "F_FUNGUS"
              found_result <- found
              found <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
            }
            # (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome) ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome)\n")
            }
            x_strip <- a.x_backup %>% strsplit("[ .]") %>% unlist()
            if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
              for (i in 2:(length(x_strip))) {
                x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- lookup(mo == found)
                  # uncertainty level 2 only if searched part contains a space (otherwise it will be found with lvl 3)
                  if (x_strip_collapsed %like_case% " ") {
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result))
                    return(found)
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
            x_strip <- a.x_backup %>% strsplit("[ .]") %>% unlist()
            if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
              for (i in 2:(length(x_strip))) {
                x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- lookup(mo == found)
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result))
                  return(found)
                }
              }
            }
            # (11) try to strip off one element from end and check the remains (any text size) ----
            # (this is in fact 7 but without nchar limit of >=6)
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (11) try to strip off one element from end and check the remains (any text size)\n")
            }
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                x_strip_collapsed <- paste(x_strip[seq_len(length(x_strip) - i)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- lookup(mo == found)
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result))
                  return(found)
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
            found <- lookup(fullname_lower %like_case% f.x_withspaces_end_only, column = "mo")
            if (!is.na(found) & nchar(g.x_backup_without_spp) >= 6) {
              found_result <- lookup(mo == found)
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result))
              return(found)
            }
          }
          
          
          # didn't found in uncertain results too
          return(NA_character_)
        }
        
        # uncertain results
        x[i] <- uncertain_fn(a.x_backup = a.x_backup, 
                             b.x_trimmed = b.x_trimmed,
                             d.x_withspaces_start_end = d.x_withspaces_start_end,
                             e.x_withspaces_start_only = e.x_withspaces_start_only, 
                             f.x_withspaces_end_only = f.x_withspaces_end_only,
                             g.x_backup_without_spp = g.x_backup_without_spp,
                             uncertain.reference_data_to_use = MO_lookup[which(MO_lookup$prevalence %in% c(1, 2)), ])
        if (!empty_result(x[i])) {
          return(x[i])
        }
        x[i] <- uncertain_fn(a.x_backup = a.x_backup, 
                             b.x_trimmed = b.x_trimmed,
                             d.x_withspaces_start_end = d.x_withspaces_start_end,
                             e.x_withspaces_start_only = e.x_withspaces_start_only, 
                             f.x_withspaces_end_only = f.x_withspaces_end_only,
                             g.x_backup_without_spp = g.x_backup_without_spp,
                             uncertain.reference_data_to_use = MO_lookup[which(MO_lookup$prevalence == 3), ])
        if (!empty_result(x[i])) {
          return(x[i])
        }
        
        # didn't found any
        return(NA_character_)
      }
      
      # CHECK ALL IN ONE GO ----
      x[i] <- check_per_prevalence(data_to_check = MO_lookup,
                                   data.old_to_check = MO.old_lookup,
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp = x_backup_without_spp[i],
                                   h.x_species = x_species[i],
                                   i.x_trimmed_species = x_trimmed_species[i])
      if (!empty_result(x[i])) {
        next
      }
      
      
      # no results found: make them UNKNOWN ----
      x[i] <- lookup(mo == "UNKNOWN")
      if (initial_search == TRUE) {
        failures <- c(failures, x_backup[i])
      }
    }
    
    if (initial_search == TRUE) {
      close(progress)
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
    x_input_clean <- trimws2(x_input)
    total_failures <- length(x_input_clean[as.character(x_input_clean) %in% as.character(failures) & !x_input %in% c(NA, NULL, NaN)])
    total_n <- length(x_input[!x_input %in% c(NA, NULL, NaN)])
    msg <- paste0(nr2char(n_distinct(failures)), " unique ", plural[1],
                  " (covering ", percentage(total_failures / total_n),
                  ") could not be coerced and ", plural[3], " considered 'unknown'")
    if (n_distinct(failures) <= 10) {
      msg <- paste0(msg, ": ", paste('"', unique(failures), '"', sep = "", collapse = ", "))
    }
    msg <- paste0(msg,  ".\nUse mo_failures() to review ", plural[2], ". Edit the `allow_uncertain` parameter if needed (see ?as.mo).")
    warning(font_red(paste0("\n", msg)),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  # handling uncertainties ----
  if (NROW(uncertainties) > 0 & initial_search == TRUE) {
    options(mo_uncertainties = as.list(distinct(uncertainties, input, .keep_all = TRUE)))
    
    plural <- c("", "it", "was")
    if (NROW(uncertainties) > 1) {
      plural <- c("s", "them", "were")
    }
    msg <- paste0("Result", plural[1], " of ", nr2char(NROW(uncertainties)), " value", plural[1],
                  " ", plural[3], " guessed with uncertainty. Use mo_uncertainties() to review ", plural[2], ".")
    warning(font_red(paste0("\n", msg)),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  
  # Becker ----
  if (Becker == TRUE | Becker == "all") {
    # See Source. It's this figure:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4187637/figure/F3/
    MOs_staph <- MO_lookup[which(MO_lookup$genus == "Staphylococcus"), ]
    CoNS <- MOs_staph[which(MOs_staph$species %in% c("arlettae", "auricularis", "capitis",
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
                            | (MOs_staph$species == "schleiferi" & MOs_staph$subspecies %in% c("schleiferi", ""))),
                      property]
    CoPS <- MOs_staph[which(MOs_staph$species %in% c("simiae", "agnetis",
                                                     "delphini", "lutrae",
                                                     "hyicus", "intermedius",
                                                     "pseudintermedius", "pseudointermedius",
                                                     "schweitzeri", "argenteus")
                            | (MOs_staph$species == "schleiferi" & MOs_staph$subspecies == "coagulans")),
                      property]
    
    # warn when species found that are not in Becker (2014, PMID 25278577) and Becker (2019, PMID 30872103)
    post_Becker <- c("argensis", "caeli", "cornubiensis", "edaphicus")
    if (any(x %in% MOs_staph[which(MOs_staph$species %in% post_Becker), property])) {
      
      warning("Becker ", font_italic("et al."), " (2014, 2019) does not contain these species named after their publication: ",
              font_italic(paste("S.",
                                sort(mo_species(unique(x[x %in% MOs_staph[which(MOs_staph$species %in% post_Becker), property]]))),
                                collapse = ", ")),
              ".",
              call. = FALSE,
              immediate. = TRUE)
    }
    
    x[x %in% CoNS] <- lookup(mo == "B_STPHY_CONS")
    x[x %in% CoPS] <- lookup(mo == "B_STPHY_COPS")
    if (Becker == "all") {
      x[x %in% lookup(fullname %like_case% "^Staphylococcus aureus", n = Inf)] <- lookup(mo == "B_STPHY_COPS")
    }
  }
  
  # Lancefield ----
  if (Lancefield == TRUE | Lancefield == "all") {
    # group A - S. pyogenes
    x[x %in% lookup(genus == "Streptococcus" & species == "pyogenes", n = Inf)] <- lookup(fullname == "Streptococcus group A")
    # group B - S. agalactiae
    x[x %in% lookup(genus == "Streptococcus" & species == "agalactiae", n = Inf)] <- lookup(fullname == "Streptococcus group B")
    # group C
    x[x %in% lookup(genus == "Streptococcus" &
                      species %in% c("equisimilis", "equi", "zooepidemicus", "dysgalactiae"),
                    n = Inf)] <- lookup(fullname == "Streptococcus group C")
    if (Lancefield == "all") {
      # all Enterococci
      x[x %in% lookup(genus == "Enterococcus", n = Inf)] <- lookup(fullname == "Streptococcus group D")
    }
    # group F - S. anginosus
    x[x %in% lookup(genus == "Streptococcus" & species == "anginosus", n = Inf)] <- lookup(fullname == "Streptococcus group F")
    # group H - S. sanguinis
    x[x %in% lookup(genus == "Streptococcus" & species == "sanguinis", n = Inf)] <- lookup(fullname == "Streptococcus group H")
    # group K - S. salivarius
    x[x %in% lookup(genus == "Streptococcus" & species == "salivarius", n = Inf)] <- lookup(fullname == "Streptococcus group K")
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
  
  # super fast using base::match() which is a lot faster than base::merge()
  x <- df_found$found[match(df_input$input, df_found$input)]
  
  if (property == "mo") {
    x <- to_class_mo(x)
  }

  if (length(mo_renamed()) > 0) {
    print(mo_renamed())
  }
  
  if (old_mo_warning == TRUE & property != "mo") {
    warning("The input contained old microorganism IDs from previous versions of this package.\nPlease use `as.mo()` on these old IDs to transform them to the new format.\nSUPPORT FOR THIS WILL BE DROPPED IN A FUTURE VERSION.", call. = FALSE)
  }
  
  x
}

empty_result <- function(x) {
  all(x %in% c(NA, "UNKNOWN"))
}

was_renamed <- function(name_old, name_new, ref_old = "", ref_new = "", mo = "") {
  newly_set <- data.frame(old_name = name_old, 
                          old_ref = ref_old,
                          new_name = name_new,
                          new_ref = ref_new,
                          mo = mo, 
                          stringsAsFactors = FALSE)
  already_set <- getOption("mo_renamed")
  if (!is.null(already_set)) {
    options(mo_renamed = rbind(already_set, newly_set))
  } else {
    options(mo_renamed = newly_set)
  }
}

format_uncertainty_as_df <- function(uncertainty_level,
                                     input,
                                     result_mo) {
  if (!is.null(getOption("mo_renamed_last_run", default = NULL))) {
    # was found as a renamed mo
    df <- data.frame(uncertainty = uncertainty_level,
                     input = input,
                     fullname = getOption("mo_renamed_last_run"),
                     renamed_to = MO_lookup[which(MO_lookup$mo == result_mo), "fullname"][1],
                     mo = result_mo,
                     stringsAsFactors = FALSE)
    options(mo_renamed_last_run = NULL)
  } else {
    df <- data.frame(uncertainty = uncertainty_level,
                     input = input,
                     fullname = MO_lookup[which(MO_lookup$mo == result_mo), "fullname"][1],
                     renamed_to = NA_character_,
                     mo = result_mo,
                     stringsAsFactors = FALSE)
  }
  df
}

#' @exportMethod print.mo
#' @export
#' @noRd
print.mo <- function(x, ...) {
  cat("Class <mo>\n")
  x_names <- names(x)
  x <- as.character(x)
  names(x) <- x_names
  print.default(x, quote = FALSE)
}

#' @exportMethod summary.mo
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- as.mo(object) # force again, could be mo from older pkg version
  top <- as.data.frame(table(x), responseName = "n", stringsAsFactors = FALSE)
  top_3 <- top[order(-top$n), 1][1:3]
  c("Class" = "mo",
    "<NA>" = length(x[is.na(x)]),
    "Unique" = n_distinct(x[!is.na(x)]),
    "#1" = top_3[1],
    "#2" = top_3[2],
    "#3" = top_3[3])
}

#' @exportMethod as.data.frame.mo
#' @export
#' @noRd
as.data.frame.mo <- function(x, ...) {
  nm <- deparse1(substitute(x))
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(as.mo(x), ..., nm = nm)
  } else {
    as.data.frame.vector(as.mo(x), ...)
  }
}

#' @exportMethod [.mo
#' @export
#' @noRd
"[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [[.mo
#' @export
#' @noRd
"[[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [<-.mo
#' @export
#' @noRd
"[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  # must only contain valid MOs
  class_integrity_check(y, "microorganism code", c(as.character(microorganisms$mo), 
                                                   as.character(microorganisms.translation$mo_old)))
}
#' @exportMethod [[<-.mo
#' @export
#' @noRd
"[[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  # must only contain valid MOs
  class_integrity_check(y, "microorganism code", c(as.character(microorganisms$mo), 
                                                   as.character(microorganisms.translation$mo_old)))
}
#' @exportMethod c.mo
#' @export
#' @noRd
c.mo <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  # must only contain valid MOs
  class_integrity_check(y, "microorganism code", c(as.character(microorganisms$mo), 
                                                   as.character(microorganisms.translation$mo_old)))
}

#' @rdname as.mo
#' @export
mo_failures <- function() {
  getOption("mo_failures")
}

#' @rdname as.mo
#' @export
mo_uncertainties <- function() {
  if (is.null(getOption("mo_uncertainties"))) {
    return(NULL)
  }
  structure(.Data = as.data.frame(getOption("mo_uncertainties"), stringsAsFactors = FALSE),
            class = c("mo_uncertainties", "data.frame"))
}

#' @exportMethod print.mo_uncertainties
#' @export
#' @noRd
print.mo_uncertainties <- function(x, ...) {
  if (NROW(x) == 0) {
    return(NULL)
  }
  cat(paste0(font_bold(nr2char(nrow(x)), paste0("unique result", ifelse(nrow(x) > 1, "s", ""), " guessed with uncertainty:")),
             "\n(1 = ", font_green("renamed/misspelled"),
             ", 2 = ", font_yellow("uncertain"),
             ", 3 = ", font_red("very uncertain"), ")\n"))
  
  msg <- ""
  for (i in seq_len(nrow(x))) {
    if (x[i, "uncertainty"] == 1) {
      colour1 <- font_green
      colour2 <- function(...) font_green_bg(font_white(...))
    } else if (x[i, "uncertainty"] == 2) {
      colour1 <- font_yellow
      colour2 <- function(...) font_yellow_bg(font_black(...))
    } else {
      colour1 <- font_red
      colour2 <- function(...) font_red_bg(font_white(...))
    }
    msg <- paste(msg,
                 paste0(colour2(paste0(" [", x[i, "uncertainty"], "] ")), ' "', x[i, "input"], '" -> ',
                        colour1(paste0(font_italic(x[i, "fullname"]),
                                       ifelse(!is.na(x[i, "renamed_to"]), paste(", renamed to", font_italic(x[i, "renamed_to"])), ""),
                                       " (", x[i, "mo"],
                                       ", score: ", percentage(levenshtein_fraction(x[i, "input"], x[i, "fullname"]), digits = 1),
                                       ")"))),
                 sep = "\n")
  }
  cat(msg)
}

#' @rdname as.mo
#' @export
mo_renamed <- function() {
  items <- getOption("mo_renamed")
  if (is.null(items)) {
    items <- data.frame()
  } else {
    items <- distinct(items, old_name, .keep_all = TRUE)
  }
  structure(.Data = items,
            class = c("mo_renamed", "data.frame"))
}

#' @exportMethod print.mo_renamed
#' @export
#' @noRd
print.mo_renamed <- function(x, ...) {
  if (NROW(x) == 0) {
    return(invisible())
  }
  for (i in seq_len(nrow(x))) {
    message(font_blue(paste0("NOTE: ", 
                             font_italic(x$old_name[i]), ifelse(x$old_ref[i] %in% c("", NA), "", 
                                                                paste0(" (",  gsub("et al.", font_italic("et al."), x$old_ref[i]), ")")),
                             " was renamed ", 
                             ifelse(as.integer(gsub("[^0-9]", "", x$new_ref[i])) < as.integer(gsub("[^0-9]", "", x$old_ref[i])),
                                    font_bold("back to "),
                                    ""),
                             font_italic(x$new_name[i]), ifelse(x$new_ref[i] %in% c("", NA), "", 
                                                                paste0(" (",  gsub("et al.", font_italic("et al."), x$new_ref[i]), ")")),
                             " [", x$mo[i], "]")))
  }
}

nr2char <- function(x) {
  if (x %in% c(1:10)) {
    v <- c("one" = 1, "two" = 2, "three" = 3, "four" = 4, "five" = 5,
           "six" = 6, "seven" = 7, "eight" = 8, "nine" = 9, "ten" = 10)
    names(v[x])
  } else {
    x
  }
}

unregex <- function(x) {
  gsub("[^a-zA-Z0-9 -]", "", x)
}

translate_allow_uncertain <- function(allow_uncertain) {
  if (isTRUE(allow_uncertain)) {
    # default to uncertainty level 2
    allow_uncertain <- 2
  } else {
    allow_uncertain[tolower(allow_uncertain) == "none"] <- 0
    allow_uncertain[tolower(allow_uncertain) == "all"] <- 3
    allow_uncertain <- as.integer(allow_uncertain)
    if (!allow_uncertain %in% c(0:3)) {
      stop('`allow_uncertain` must be a number between 0 (or "none") and 3 (or "all"), or TRUE (= 2) or FALSE (= 0).', call. = FALSE)
    }
  }
  allow_uncertain
}

get_mo_failures_uncertainties_renamed <- function() {
  remember <- list(failures = getOption("mo_failures"),
                   uncertainties = getOption("mo_uncertainties"),
                   renamed = getOption("mo_renamed"))
  # empty them, otherwise mo_shortname("Chlamydophila psittaci") will give 3 notes
  options("mo_failures" = NULL)
  options("mo_uncertainties" = NULL)
  options("mo_renamed" = NULL)
  remember
}

load_mo_failures_uncertainties_renamed <- function(metadata) {
  options("mo_failures" = metadata$failures)
  options("mo_uncertainties" = metadata$uncertainties)
  options("mo_renamed" = metadata$renamed)
}

levenshtein_fraction <- function(input, output) {
  levenshtein <- double(length = length(input))
  for (i in seq_len(length(input))) {
    # determine Levenshtein distance, but maximise to nchar of output
    levenshtein[i] <- base::min(base::as.double(utils::adist(input[i], output[i], ignore.case = TRUE)),
                                base::nchar(output[i]))
  }
  # self-made score between 0 and 1 (for % certainty, so 0 means huge distance, 1 means no distance)
  (base::nchar(output) - 0.5 * levenshtein) / nchar(output)
}

trimws2 <- function(x) {
  trimws(gsub("[\\s]+", " ", x, perl = TRUE))
}

parse_and_convert <- function(x) {
  tryCatch({
    if (!is.null(dim(x))) {
      if (NCOL(x) > 2) {
        stop("A maximum of two columns is allowed.", call. = FALSE)
      } else if (NCOL(x) == 2) {
        # support tidyverse selection like: df %>% select(colA, colB)
        # paste these columns together
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        colnames(x) <- c("A", "B")
        x <- paste(x$A, x$B)
      } else {
        # support tidyverse selection like: df %>% select(colA)
        x <- as.data.frame(x, stringsAsFactors = FALSE)[[1]]
      }
    }
    x[is.null(x)] <- NA
    parsed <- iconv(x, to = "UTF-8")
    parsed[is.na(parsed) & !is.na(x)] <- iconv(x[is.na(parsed) & !is.na(x)], from = "Latin1", to = "ASCII//TRANSLIT")
    parsed <- gsub('"', "", parsed, fixed = TRUE)
  }, error = function(e) stop(e$message, call. = FALSE)) # this will also be thrown when running `as.mo(no_existing_object)`
  parsed
}

left_join_MO_lookup <- function(x, ...) {
  left_join(x = x, y = MO_lookup, ...)
}
left_join_MO.old_lookup <- function(x, ...) {
  left_join(x = x, y = MO.old_lookup, ...)
}

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

#' Transform Arbitrary Input to Valid Microbial Taxonomy
#'
#' Use this function to get a valid microorganism code ([`mo`]) based on arbitrary user input. Determination is done using intelligent rules and the complete taxonomic tree of the kingdoms `r vector_and(unique(microorganisms$kingdom[which(!grepl("(unknown|Fungi)", microorganisms$kingdom))]), quotes = FALSE)`, and most microbial species from the kingdom Fungi (see *Source*). The input can be almost anything: a full name (like `"Staphylococcus aureus"`), an abbreviated name (such as `"S. aureus"`), an abbreviation known in the field (such as `"MRSA"`), or just a genus. See *Examples*.
#' @param x A [character] vector or a [data.frame] with one or two columns.
#' @param Becker A [logical] to indicate whether staphylococci should be categorised into coagulase-negative staphylococci ("CoNS") and coagulase-positive staphylococci ("CoPS") instead of their own species, according to Karsten Becker *et al.* (see *Source*). Please see *Details* for a full list of staphylococcal species that will be converted.
#'
#' This excludes *Staphylococcus aureus* at default, use `Becker = "all"` to also categorise *S. aureus* as "CoPS".
#' @param Lancefield A [logical] to indicate whether a beta-haemolytic *Streptococcus* should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield (see *Source*). These streptococci will be categorised in their first group, e.g. *Streptococcus dysgalactiae* will be group C, although officially it was also categorised into groups G and L. . Please see *Details* for a full list of streptococcal species that will be converted.
#'
#' This excludes enterococci at default (who are in group D), use `Lancefield = "all"` to also categorise all enterococci as group D.
#' @param minimum_matching_score A numeric value to set as the lower limit for the [MO matching score][mo_matching_score()]. When left blank, this will be determined automatically based on the character length of `x`, its [taxonomic kingdom][microorganisms] and [human pathogenicity][mo_matching_score()].
#' @param keep_synonyms A [logical] to indicate if old, previously valid taxonomic names must be preserved and not be corrected to currently accepted names. The default is `FALSE`, which will return a note if old taxonomic names were processed. The default can be set with the package option [`AMR_keep_synonyms`][AMR-options], i.e. `options(AMR_keep_synonyms = TRUE)` or `options(AMR_keep_synonyms = FALSE)`.
#' @param reference_df A [data.frame] to be used for extra reference when translating `x` to a valid [`mo`]. See [set_mo_source()] and [get_mo_source()] to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ignore_pattern A Perl-compatible [regular expression][base::regex] (case-insensitive) of which all matches in `x` must return `NA`. This can be convenient to exclude known non-relevant input and can also be set with the package option [`AMR_ignore_pattern`][AMR-options], e.g. `options(AMR_ignore_pattern = "(not reported|contaminated flora)")`.
#' @param cleaning_regex A Perl-compatible [regular expression][base::regex] (case-insensitive) to clean the input of `x`. Every matched part in `x` will be removed. At default, this is the outcome of [mo_cleaning_regex()], which removes texts between brackets and texts such as "species" and "serovar". The default can be set with the package option [`AMR_cleaning_regex`][AMR-options].
#' @param only_fungi A [logical] to indicate if only fungi must be found, making sure that e.g. misspellings always return records from the kingdom of Fungi. This can be set globally for [all microorganism functions][mo_property()] with the package option [`AMR_only_fungi`][AMR-options], i.e. `options(AMR_only_fungi = TRUE)`.
#' @param language Language to translate text like "no growth", which defaults to the system language (see [get_AMR_locale()]).
#' @param info A [logical] to indicate that info must be printed, e.g. a progress bar when more than 25 items are to be coerced, or a list with old taxonomic names. The default is `TRUE` only in interactive mode.
#' @param ... Other arguments passed on to functions.
#' @rdname as.mo
#' @aliases mo
#' @details
#' A microorganism (MO) code from this package (class: [`mo`]) is human-readable and typically looks like these examples:
#'
#' ```
#'   Code               Full name
#'   ---------------    --------------------------------------
#'   B_KLBSL            Klebsiella
#'   B_KLBSL_PNMN       Klebsiella pneumoniae
#'   B_KLBSL_PNMN_RHNS  Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |    |
#'   |   |    |    |
#'   |   |    |    \---> subspecies, a 3-5 letter acronym
#'   |   |    \----> species, a 3-6 letter acronym
#'   |   \----> genus, a 4-8 letter acronym
#'   \----> kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                   C (Chromista), F (Fungi), PL (Plantae),
#'                   P (Protozoa)
#' ```
#'
#' Values that cannot be coerced will be considered 'unknown' and will return the MO code `UNKNOWN` with a warning.
#'
#' Use the [`mo_*`][mo_property()] functions to get properties based on the returned code, see *Examples*.
#'
#' The [as.mo()] function uses a novel and scientifically validated (\doi{10.18637/jss.v104.i03}) matching score algorithm (see *Matching Score for Microorganisms* below) to match input against the [available microbial taxonomy][microorganisms] in this package. This implicates that e.g. `"E. coli"` (a microorganism highly prevalent in humans) will return the microbial ID of *Escherichia coli* and not *Entamoeba coli* (a microorganism less prevalent in humans), although the latter would alphabetically come first.
#'
#' ### Coping with Uncertain Results
#'
#' Results of non-exact taxonomic input are based on their [matching score][mo_matching_score()]. The lowest allowed score can be set with the `minimum_matching_score` argument. At default this will be determined based on the character length of the input, the [taxonomic kingdom][microorganisms], and the [human pathogenicity][mo_matching_score()] of the taxonomic outcome. If values are matched with uncertainty, a message will be shown to suggest the user to inspect the results with [mo_uncertainties()], which returns a [data.frame] with all specifications.
#'
#' To increase the quality of matching, the `cleaning_regex` argument is used to clean the input. This must be a [regular expression][base::regex] that matches parts of the input that should be removed before the input is matched against the [available microbial taxonomy][microorganisms]. It will be matched Perl-compatible and case-insensitive. The default value of `cleaning_regex` is the outcome of the helper function [mo_cleaning_regex()].
#'
#' There are three helper functions that can be run after using the [as.mo()] function:
#' - Use [mo_uncertainties()] to get a [data.frame] that prints in a pretty format with all taxonomic names that were guessed. The output contains the matching score for all matches (see *Matching Score for Microorganisms* below).
#' - Use [mo_failures()] to get a [character] [vector] with all values that could not be coerced to a valid value.
#' - Use [mo_renamed()] to get a [data.frame] with all values that could be coerced based on old, previously accepted taxonomic names.
#'
#' ### For Mycologists
#'
#' The [matching score algorithm][mo_matching_score()] gives precedence to bacteria over fungi. If you are only analysing fungi, be sure to use `only_fungi = TRUE`, or better yet, add this to your code and run it once every session:
#'
#' ```r
#' options(AMR_only_fungi = TRUE)
#' ```
#'
#' This will make sure that no bacteria or other 'non-fungi' will be returned by [as.mo()], or any of the [`mo_*`][mo_property()] functions.
#'
#' ### Coagulase-negative and Coagulase-positive Staphylococci
#'
#' With `Becker = TRUE`, the following staphylococci will be converted to their corresponding coagulase group:
#'
#' * Coagulase-negative: `r vector_and(gsub("Staphylococcus", "S.", mo_name(MO_CONS[MO_CONS != "B_STPHY_CONS"], keep_synonyms = TRUE)), quotes = "*")`
#' * Coagulase-positive: `r vector_and(gsub("Staphylococcus", "S.", mo_name(MO_COPS[MO_COPS != "B_STPHY_COPS"], keep_synonyms = TRUE)), quotes = "*")`
#'
#' This is based on:
#'
#' * Becker K *et al.* (2014). **Coagulase-Negative Staphylococci.** *Clin Microbiol Rev.* 27(4): 870-926; \doi{10.1128/CMR.00109-13}
#' * Becker K *et al.* (2019). **Implications of identifying the recently defined members of the *S. aureus* complex, *S. argenteus* and *S. schweitzeri*: A position paper of members of the ESCMID Study Group for staphylococci and Staphylococcal Diseases (ESGS).** *Clin Microbiol Infect*; \doi{10.1016/j.cmi.2019.02.028}
#' * Becker K *et al.* (2020). **Emergence of coagulase-negative staphylococci.** *Expert Rev Anti Infect Ther.* 18(4):349-366; \doi{10.1080/14787210.2020.1730813}
#'
#' For newly named staphylococcal species, such as *S. brunensis* (2024) and *S. shinii* (2023), we looked up the scientific reference to make sure the species are considered for the correct coagulase group.
#'
#' ### Lancefield Groups in Streptococci
#'
#' With `Lancefield = TRUE`, the following streptococci will be converted to their corresponding Lancefield group:
#'
#' * `r paste(apply(aggregate(mo_name ~ mo_group_name, data = microorganisms.groups[microorganisms.groups$mo_group_name %like_case% "Streptococcus Group [A-Z]$", ], FUN = function(x) vector_and(gsub("Streptococcus", "S.", x, fixed = TRUE), quotes = "*", sort = TRUE)), 1, function(row) paste(row["mo_group_name"], ": ", row["mo_name"], sep = "")), collapse = "\n* ")`
#'
#' This is based on:
#'
#' * Lancefield RC (1933). **A serological differentiation of human and other groups of hemolytic streptococci.** *J Exp Med.* 57(4): 571-95; \doi{10.1084/jem.57.4.571}
#'
#' @inheritSection mo_matching_score Matching Score for Microorganisms
#'
#  (source as a section here, so it can be inherited by other man pages)
#' @section Source:
#' * Berends MS *et al.* (2022). **AMR: An R Package for Working with Antimicrobial Resistance Data**. *Journal of Statistical Software*, 104(3), 1-31; \doi{10.18637/jss.v104.i03}
#' * `r TAXONOMY_VERSION$LPSN$citation` Accessed from <`r TAXONOMY_VERSION$LPSN$url`> on `r documentation_date(TAXONOMY_VERSION$LPSN$accessed_date)`.
#' * `r TAXONOMY_VERSION$MycoBank$citation` Accessed from <`r TAXONOMY_VERSION$MycoBank$url`> on `r documentation_date(TAXONOMY_VERSION$MycoBank$accessed_date)`.
#' * `r TAXONOMY_VERSION$GBIF$citation` Accessed from <`r TAXONOMY_VERSION$GBIF$url`> on `r documentation_date(TAXONOMY_VERSION$GBIF$accessed_date)`.
#' * `r TAXONOMY_VERSION$BacDive$citation` Accessed from <`r TAXONOMY_VERSION$BacDive$url`> on `r documentation_date(TAXONOMY_VERSION$BacDive$accessed_date)`.
#' * `r TAXONOMY_VERSION$SNOMED$citation` URL: <`r TAXONOMY_VERSION$SNOMED$url`>
#' * Bartlett A *et al.* (2022). **A comprehensive list of bacterial pathogens infecting humans** *Microbiology* 168:001269; \doi{10.1099/mic.0.001269}
#' @export
#' @return A [character] [vector] with additional class [`mo`]
#' @seealso [microorganisms] for the [data.frame] that is being used to determine ID's.
#'
#' The [`mo_*`][mo_property()] functions (such as [mo_genus()], [mo_gramstain()]) to get properties based on the returned code.
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' \donttest{
#' # These examples all return "B_STPHY_AURS", the ID of S. aureus:
#' as.mo(c(
#'   "sau", # WHONET code
#'   "stau",
#'   "STAU",
#'   "staaur",
#'   "S. aureus",
#'   "S aureus",
#'   "Sthafilokkockus aureus", # handles incorrect spelling
#'   "Staphylococcus aureus (MRSA)",
#'   "MRSA", # Methicillin Resistant S. aureus
#'   "VISA", # Vancomycin Intermediate S. aureus
#'   "VRSA", # Vancomycin Resistant S. aureus
#'   115329001 # SNOMED CT code
#' ))
#'
#' # Dyslexia is no problem - these all work:
#' as.mo(c(
#'   "Ureaplasma urealyticum",
#'   "Ureaplasma urealyticus",
#'   "Ureaplasmium urealytica",
#'   "Ureaplazma urealitycium"
#' ))
#'
#' # input will get cleaned up with the input given in the `cleaning_regex` argument,
#' # which defaults to `mo_cleaning_regex()`:
#' cat(mo_cleaning_regex(), "\n")
#'
#' as.mo("Streptococcus group A")
#'
#' as.mo("S. epidermidis") # will remain species: B_STPHY_EPDR
#' as.mo("S. epidermidis", Becker = TRUE) # will not remain species: B_STPHY_CONS
#'
#' as.mo("S. pyogenes") # will remain species: B_STRPT_PYGN
#' as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRPA
#'
#' # All mo_* functions use as.mo() internally too (see ?mo_property):
#' mo_genus("E. coli")
#' mo_gramstain("ESCO")
#' mo_is_intrinsic_resistant("ESCCOL", ab = "vanco")
#' }
as.mo <- function(x,
                  Becker = FALSE,
                  Lancefield = FALSE,
                  minimum_matching_score = NULL,
                  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
                  reference_df = get_mo_source(),
                  ignore_pattern = getOption("AMR_ignore_pattern", NULL),
                  cleaning_regex = getOption("AMR_cleaning_regex", mo_cleaning_regex()),
                  only_fungi = getOption("AMR_only_fungi", FALSE),
                  language = get_AMR_locale(),
                  info = interactive(),
                  ...) {
  meet_criteria(x, allow_class = c("mo", "data.frame", "list", "character", "numeric", "integer", "factor"), allow_NA = TRUE)
  meet_criteria(Becker, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(Lancefield, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(minimum_matching_score, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)
  meet_criteria(reference_df, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(ignore_pattern, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(cleaning_regex, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(only_fungi, allow_class = "logical", has_length = 1)
  language <- validate_language(language)
  meet_criteria(info, allow_class = "logical", has_length = 1)

  add_MO_lookup_to_AMR_env()

  if (tryCatch(all(x %in% c(AMR_env$MO_lookup$mo, NA)), error = function(e) FALSE) &&
    isFALSE(Becker) &&
    isFALSE(Lancefield) &&
    isTRUE(keep_synonyms)) {
    # don't look into valid MO codes, just return them
    # is.mo() won't work - MO codes might change between package versions
    return(set_clean_class(x, new_class = c("mo", "character")))
  }

  # start off with replaced language-specific non-ASCII characters with ASCII characters
  x <- parse_and_convert(x)
  # replace mo codes used in older package versions
  x <- replace_old_mo_codes(x, property = "mo")
  # ignore cases that match the ignore pattern
  x <- replace_ignore_pattern(x, ignore_pattern)

  x_lower <- tolower(x)

  # WHONET: xxx = no growth
  x[x_lower %in% c("", "xxx", "na", "nan")] <- NA_character_

  out <- rep(NA_character_, length(x))

  # below we use base R's match(), known for powering '%in%', and incredibly fast!

  # From reference_df ----
  reference_df <- repair_reference_df(reference_df)
  if (!is.null(reference_df)) {
    out[x %in% reference_df[[1]]] <- reference_df[[2]][match(x[x %in% reference_df[[1]]], reference_df[[1]])]
  }
  # From MO code ----
  out[is.na(out) & toupper(x) %in% AMR_env$MO_lookup$mo] <- toupper(x[is.na(out) & toupper(x) %in% AMR_env$MO_lookup$mo])
  # From full name ----
  out[is.na(out) & x_lower %in% AMR_env$MO_lookup$fullname_lower] <- AMR_env$MO_lookup$mo[match(x_lower[is.na(out) & x_lower %in% AMR_env$MO_lookup$fullname_lower], AMR_env$MO_lookup$fullname_lower)]
  # one exception: "Fungi" matches the kingdom, but instead it should return the 'unknown' code for fungi
  out[out == "F_[KNG]_FUNGI"] <- "F_FUNGUS"
  # From known codes ----
  ind <- is.na(out) & toupper(x) %in% AMR::microorganisms.codes$code
  out[ind] <- AMR::microorganisms.codes$mo[match(toupper(x)[ind], AMR::microorganisms.codes$code)]
  if (length(which(ind)) > 0 && isTRUE(info) && message_not_thrown_before("as.mo_microorganisms.codes", is.na(out), toupper(x))) {
    message_(
      "Retrieved value", ifelse(sum(ind) > 1, "s", ""),
      " from the `microorganisms.codes` data set for ", vector_and(toupper(x)[ind]), "."
    )
  }
  # From SNOMED ----
  # based on this extremely fast gem: https://stackoverflow.com/a/11002456/4575331
  snomeds <- unlist(AMR_env$MO_lookup$snomed)
  snomeds <- snomeds[!is.na(snomeds)]
  out[is.na(out) & x %in% snomeds] <- AMR_env$MO_lookup$mo[rep(seq_along(AMR_env$MO_lookup$snomed), vapply(FUN.VALUE = double(1), AMR_env$MO_lookup$snomed, length))[match(x[is.na(out) & x %in% snomeds], snomeds)]]
  # From other familiar output ----
  # such as Salmonella groups, colloquial names, etc.
  out[is.na(out)] <- convert_colloquial_input(x[is.na(out)])
  # From previous hits in this session ----
  old <- out
  out[is.na(out) & paste(x, minimum_matching_score, only_fungi) %in% AMR_env$mo_previously_coerced$x] <- AMR_env$mo_previously_coerced$mo[match(paste(x, minimum_matching_score, only_fungi)[is.na(out) & paste(x, minimum_matching_score, only_fungi) %in% AMR_env$mo_previously_coerced$x], AMR_env$mo_previously_coerced$x)]
  new <- out
  if (isTRUE(info) && message_not_thrown_before("as.mo", old, new, entire_session = TRUE) && any(is.na(old) & !is.na(new), na.rm = TRUE)) {
    message_(
      "Returning previously coerced value", ifelse(sum(is.na(old) & !is.na(new)) > 1, "s", ""),
      " for ", vector_and(x[is.na(old) & !is.na(new)]), ". Run `mo_reset_session()` to reset this. This note will be shown once per session for this input."
    )
  }

  # For all other input ----
  if (any(is.na(out) & !is.na(x))) {
    # reset uncertainties
    AMR_env$mo_uncertainties <- AMR_env$mo_uncertainties[0, ]
    AMR_env$mo_failures <- NULL

    # Laboratory systems: remove (translated) entries like "no growth", "not E. coli", etc.
    x[trimws2(x) %like% translate_AMR("no .*growth", language = language)] <- NA_character_
    x[trimws2(x) %like% paste0("^(", translate_AMR("no|not", language = language), ") ")] <- NA_character_

    # groups are in our taxonomic table with a capital G
    x <- gsub(" group( |$)", " Group\\1", x, perl = TRUE)

    # convert translations
    x[x %like_case% "enter[o\u00F6]?[ck]o[ck](ken)?$"] <- gsub("(.* )?enter[o\u00F6]?[ck]o[ck](ken)?$", "enterococcus", x[x %like_case% "enter[o\u00F6]?[ck]o[ck](ken)?$"], perl = TRUE)
    x[x %like_case% "strept[o\u00F6]?[ck]o[ck](ken)?$"] <- gsub("(.* )?strept[o\u00F6]?[ck]o[ck](ken)?$", "streptococcus", x[x %like_case% "strept[o\u00F6]?[ck]o[ck](ken)?$"], perl = TRUE)
    x[x %like_case% "staph[yij]?[lo]*[ck]o[ck](ken)?$"] <- gsub("(.* )?staph[yij]?[lo]*[ck]o[ck](ken)?$", "staphylococcus", x[x %like_case% "staph[yij]?[lo]*[ck]o[ck](ken)?$"], perl = TRUE)

    # run over all unique leftovers
    x_unique <- unique(x[is.na(out) & !is.na(x)])

    # set up progress bar
    progress <- progress_ticker(n = length(x_unique), n_min = 10, print = info, title = "Converting microorganism input")
    on.exit(close(progress))

    msg <- character(0)

    MO_lookup_current <- AMR_env$MO_lookup
    if (isTRUE(only_fungi)) {
      MO_lookup_current <- MO_lookup_current[MO_lookup_current$kingdom == "Fungi", , drop = FALSE]
    }

    # run it
    x_coerced <- vapply(FUN.VALUE = character(1), x_unique, function(x_search) {
      progress$tick()

      # some required cleaning steps
      x_out <- trimws2(x_search)
      # this applies the `cleaning_regex` argument, which defaults to mo_cleaning_regex()
      x_out <- gsub(cleaning_regex, " ", x_out, ignore.case = TRUE, perl = TRUE)
      x_out <- trimws2(gsub(" +", " ", x_out, perl = TRUE))
      x_search_cleaned <- x_out
      x_out <- tolower(x_out)
      # when x_search_cleaned are only capitals (such as in codes), make them lowercase to increase matching score
      x_search_cleaned[x_search_cleaned == toupper(x_search_cleaned)] <- x_out[x_search_cleaned == toupper(x_search_cleaned)]

      # first check if cleaning led to an exact result, case-insensitive
      if (x_out %in% MO_lookup_current$fullname_lower) {
        return(as.character(MO_lookup_current$mo[match(x_out, MO_lookup_current$fullname_lower)]))
      }

      # input must not be too short
      if (nchar(x_out) < 3) {
        return("UNKNOWN")
      }

      # take out the parts, split by space
      x_parts <- strsplit(gsub("-", " ", x_out, fixed = TRUE), " ", fixed = TRUE)[[1]]
      # do a pre-match on first character (and if it contains a space, first chars of first two terms)
      if (length(x_parts) %in% c(2, 3)) {
        # for genus + species + subspecies
        if (paste(x_parts[1:2], collapse = " ") %in% MO_lookup_current$fullname_lower) {
          filtr <- which(MO_lookup_current$fullname_lower %like% paste(x_parts[1:2], collapse = " "))
        } else if (x_parts[1] %in% MO_lookup_current$genus_lower && !paste(x_parts[1:2], collapse = " ") %in% MO_lookup_current$fullname_lower) {
          # for a known genus, but unknown (sub)species
          filtr <- which(MO_lookup_current$genus_lower == x_parts[1])
          minimum_matching_score <- 0.05
        } else if (nchar(gsub("[^a-z]", "", x_parts[1], perl = TRUE)) <= 3) {
          filtr <- which(MO_lookup_current$full_first == substr(x_parts[1], 1, 1) &
            (MO_lookup_current$species_first == substr(x_parts[2], 1, 1) |
              MO_lookup_current$subspecies_first == substr(x_parts[2], 1, 1) |
              MO_lookup_current$subspecies_first == substr(x_parts[3], 1, 1)))
        } else {
          filtr <- which(MO_lookup_current$full_first == substr(x_parts[1], 1, 1) |
            MO_lookup_current$species_first == substr(x_parts[2], 1, 1) |
            MO_lookup_current$subspecies_first == substr(x_parts[2], 1, 1) |
            MO_lookup_current$subspecies_first == substr(x_parts[3], 1, 1))
        }
      } else if (length(x_parts) > 3) {
        first_chars <- paste0("(^| )[", paste(substr(x_parts, 1, 1), collapse = ""), "]")
        filtr <- which(MO_lookup_current$full_first %like_case% first_chars)
      } else if (nchar(x_out) == 3) {
        # no space and 3 characters - probably a code such as SAU or ECO
        msg <<- c(msg, paste0("Input \"", x_search, "\" was assumed to be a microorganism code - tried to match on \"", totitle(substr(x_out, 1, 1)), AMR_env$ellipsis_icon, " ", substr(x_out, 2, 3), AMR_env$ellipsis_icon, "\""))
        filtr <- which(MO_lookup_current$fullname_lower %like_case% paste0("(^| )", substr(x_out, 1, 1), ".* ", substr(x_out, 2, 3)))
      } else if (nchar(x_out) == 4) {
        # no space and 4 characters - probably a code such as STAU or ESCO
        msg <<- c(msg, paste0("Input \"", x_search, "\" was assumed to be a microorganism code - tried to match on \"", totitle(substr(x_out, 1, 2)), AMR_env$ellipsis_icon, " ", substr(x_out, 3, 4), AMR_env$ellipsis_icon, "\""))
        filtr <- which(MO_lookup_current$fullname_lower %like_case% paste0("(^| )", substr(x_out, 1, 2), ".* ", substr(x_out, 3, 4)))
      } else if (nchar(x_out) <= 6) {
        # no space and 5-6 characters - probably a code such as STAAUR or ESCCOL
        first_part <- paste0(substr(x_out, 1, 2), "[a-z]*", substr(x_out, 3, 3))
        second_part <- substr(x_out, 4, nchar(x_out))
        msg <<- c(msg, paste0("Input \"", x_search, "\" was assumed to be a microorganism code - tried to match on \"", gsub("[a-z]*", AMR_env$ellipsis_icon, totitle(first_part), fixed = TRUE), " ", second_part, AMR_env$ellipsis_icon, "\""))
        filtr <- which(MO_lookup_current$fullname_lower %like_case% paste0("(^| )", first_part, ".* ", second_part))
      } else {
        # for genus or species or subspecies
        filtr <- which(MO_lookup_current$full_first == substr(x_parts, 1, 1) |
          MO_lookup_current$species_first == substr(x_parts, 1, 1) |
          MO_lookup_current$subspecies_first == substr(x_parts, 1, 1))
      }
      if (length(filtr) == 0) {
        mo_to_search <- MO_lookup_current$fullname
      } else {
        mo_to_search <- MO_lookup_current$fullname[filtr]
      }

      AMR_env$mo_to_search <- mo_to_search
      # determine the matching score on the original search value
      m <- mo_matching_score(x = x_search_cleaned, n = mo_to_search)
      if (is.null(minimum_matching_score)) {
        minimum_matching_score_current <- min(0.6, min(10, nchar(x_search_cleaned)) * 0.08)
        # correct back for prevalence
        minimum_matching_score_current <- minimum_matching_score_current / MO_lookup_current$prevalence[match(mo_to_search, MO_lookup_current$fullname)]
        # correct back for kingdom
        minimum_matching_score_current <- minimum_matching_score_current / MO_lookup_current$kingdom_index[match(mo_to_search, MO_lookup_current$fullname)]
        minimum_matching_score_current <- pmax(minimum_matching_score_current, m)
        if (length(x_parts) > 1 && all(m <= 0.55, na.rm = TRUE)) {
          # if the highest score is 0.5, we have nothing serious - 0.5 is the lowest for pathogenic group 1
          # make everything NA so the results will get removed below
          # (we added length(x_parts) > 1 to exclude microbial codes from this rule, such as "STAU")
          m[seq_len(length(m))] <- NA_real_
        }
      } else {
        # minimum_matching_score was set, so remove everything below it
        m[m < minimum_matching_score] <- NA_real_
        minimum_matching_score_current <- minimum_matching_score
      }

      top_hits <- mo_to_search[order(m, decreasing = TRUE, na.last = NA)] # na.last = NA will remove the NAs
      if (length(top_hits) == 0) {
        warning_("No hits found for \"", x_search, "\" with minimum_matching_score = ", ifelse(is.null(minimum_matching_score), paste0("NULL (=", round(min(minimum_matching_score_current, na.rm = TRUE), 3), ")"), minimum_matching_score), ". Try setting this value lower or even to 0.", call = FALSE)
        result_mo <- NA_character_
      } else {
        result_mo <- MO_lookup_current$mo[match(top_hits[1], MO_lookup_current$fullname)]
        AMR_env$mo_uncertainties <- rbind_AMR(
          AMR_env$mo_uncertainties,
          data.frame(
            original_input = x_search,
            input = x_search_cleaned,
            fullname = top_hits[1],
            mo = result_mo,
            candidates = ifelse(length(top_hits) > 1, paste(top_hits[2:min(99, length(top_hits))], collapse = ", "), ""),
            minimum_matching_score = ifelse(is.null(minimum_matching_score), "NULL", minimum_matching_score),
            keep_synonyms = keep_synonyms,
            stringsAsFactors = FALSE
          )
        )
        # save to package env to save time for next time
        AMR_env$mo_previously_coerced <- unique(rbind_AMR(
          AMR_env$mo_previously_coerced,
          data.frame(
            x = paste(x_search, minimum_matching_score, only_fungi),
            mo = result_mo,
            stringsAsFactors = FALSE
          )
        ))
      }
      # the actual result:
      as.character(result_mo)
    })

    # remove progress bar from console
    close(progress)
    # expand from unique again
    out[is.na(out)] <- x_coerced[match(x[is.na(out)], x_unique)]

    # Throw note about uncertainties ----
    if (isTRUE(info) && NROW(AMR_env$mo_uncertainties) > 0) {
      if (message_not_thrown_before("as.mo", "uncertainties", AMR_env$mo_uncertainties$original_input)) {
        plural <- c("", "this")
        if (length(AMR_env$mo_uncertainties$original_input) > 1) {
          plural <- c("s", "these uncertainties")
        }
        if (length(AMR_env$mo_uncertainties$original_input) <= 3) {
          examples <- vector_and(
            paste0(
              '"', AMR_env$mo_uncertainties$original_input,
              '" (assumed ', italicise(AMR_env$mo_uncertainties$fullname), ")"
            ),
            quotes = FALSE
          )
        } else {
          examples <- paste0(nr2char(length(AMR_env$mo_uncertainties$original_input)), " microorganism", plural[1])
        }
        msg <- c(msg, paste0(
          "Microorganism translation was uncertain for ", examples,
          ". Run `mo_uncertainties()` to review ", plural[2], ", or use `add_custom_microorganisms()` to add custom entries."
        ))

        for (m in msg) {
          message_(m)
        }
      }
    }
  } # end of loop over all yet unknowns

  # Keep or replace synonyms ----
  out_current <- synonym_mo_to_accepted_mo(out, fill_in_accepted = FALSE)
  AMR_env$mo_renamed <- list(old = out[!is.na(out_current)])
  if (isFALSE(keep_synonyms)) {
    out[!is.na(out_current)] <- out_current[!is.na(out_current)]
    if (isTRUE(info) && length(AMR_env$mo_renamed$old) > 0) {
      print(mo_renamed(), extra_txt = " (use `keep_synonyms = TRUE` to leave uncorrected)")
    }
  } else if (is.null(getOption("AMR_keep_synonyms")) && length(AMR_env$mo_renamed$old) > 0 && message_not_thrown_before("as.mo", "keep_synonyms_warning", entire_session = TRUE)) {
    # keep synonyms is TRUE, so check if any do have synonyms
    warning_("Function `as.mo()` returned ", nr2char(length(unique(AMR_env$mo_renamed$old))), " old taxonomic name", ifelse(length(unique(AMR_env$mo_renamed$old)) > 1, "s", ""), ". Use `as.mo(..., keep_synonyms = FALSE)` to clean the input to currently accepted taxonomic names, or set the R option `AMR_keep_synonyms` to `FALSE`. This warning will be shown once per session.", call = FALSE)
  }

  # Apply Becker ----
  if (!isTRUE(only_fungi) && (isTRUE(Becker) || Becker == "all")) {
    # warn when species found that are not in:
    # - Becker et al. 2014, PMID 25278577
    # - Becker et al. 2019, PMID 30872103
    # - Becker et al. 2020, PMID 32056452

    # comment below code if all staphylococcal species are categorised as CoNS/CoPS
    post_Becker <- paste(
      "Staphylococcus",
      c("caledonicus", "canis", "durrellii", "lloydii", "ratti", "roterodami", "singaporensis", "taiwanensis")
    )
    if (any(out %in% AMR_env$MO_lookup$mo[match(post_Becker, AMR_env$MO_lookup$fullname)])) {
      if (message_not_thrown_before("as.mo", "becker")) {
        warning_("in `as.mo()`: Becker ", font_italic("et al."), " (2014, 2019, 2020) does not contain these species named after their publication: ",
          vector_and(font_italic(gsub("Staphylococcus", "S.", post_Becker, fixed = TRUE), collapse = NULL), quotes = FALSE),
          ". Categorisation to CoNS/CoPS was taken from the original scientific publication(s).",
          immediate = TRUE, call = FALSE
        )
      }
    }

    # 'MO_CONS' and 'MO_COPS' are 'mo' vectors created in R/_pre_commit_checks.R
    out[out %in% MO_CONS] <- "B_STPHY_CONS"
    out[out %in% MO_COPS] <- "B_STPHY_COPS"
    if (Becker == "all") {
      out[out == "B_STPHY_AURS"] <- "B_STPHY_COPS"
    }
  }

  # Apply Lancefield ----
  if (!isTRUE(only_fungi) && (isTRUE(Lancefield) || Lancefield == "all")) {
    # (using `%like_case%` to also match subspecies)

    # group A - S. pyogenes
    out[out %like_case% "^B_STRPT_PYGN(_|$)"] <- "B_STRPT_GRPA"
    # group B - S. agalactiae
    out[out %like_case% "^B_STRPT_AGLC(_|$)"] <- "B_STRPT_GRPB"
    # group C - all subspecies within S. dysgalactiae and S. equi (such as S. equi zooepidemicus)
    out[out %like_case% "^B_STRPT_(DYSG|EQUI)(_|$)"] <- "B_STRPT_GRPC"
    if (Lancefield == "all") {
      # group D - all enterococci
      out[out %like_case% "^B_ENTRC(_|$)"] <- "B_STRPT_GRPD"
    }
    # group F - Milleri group == S. anginosus group, which incl. S. anginosus, S. constellatus, S. intermedius
    out[out %like_case% "^B_STRPT_(ANGN|CNST|INTR)(_|$)"] <- "B_STRPT_GRPF"
    # group G - S. dysgalactiae and S. canis (though dysgalactiae is also group C and will be matched there)
    out[out %like_case% "^B_STRPT_(DYSG|CANS)(_|$)"] <- "B_STRPT_GRPG"
    # group H - S. sanguinis
    out[out %like_case% "^B_STRPT_SNGN(_|$)"] <- "B_STRPT_GRPH"
    # group K - S. salivarius, incl. S. salivarius salivarius and S. salivarius thermophilus
    out[out %like_case% "^B_STRPT_SLVR(_|$)"] <- "B_STRPT_GRPK"
    # group L - only S. dysgalactiae which is also group C & G, so ignore it here
  }

  # All unknowns ----
  out[is.na(out) & !is.na(x)] <- "UNKNOWN"
  AMR_env$mo_failures <- unique(x[out == "UNKNOWN" & !toupper(x) %in% c("UNKNOWN", "CON", "UNK") & !x %like_case% "^[(]unknown [a-z]+[)]$" & !is.na(x)])
  if (length(AMR_env$mo_failures) > 0) {
    warning_("The following input could not be coerced and was returned as \"UNKNOWN\": ", vector_and(AMR_env$mo_failures, quotes = TRUE), ".\nYou can retrieve this list with `mo_failures()`.", call = FALSE)
  }

  # Return class ----
  set_clean_class(out,
    new_class = c("mo", "character")
  )
}

# OTHER DOCUMENTED FUNCTIONS ----------------------------------------------

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  inherits(x, "mo")
}

#' @rdname as.mo
#' @export
mo_uncertainties <- function() {
  set_clean_class(AMR_env$mo_uncertainties, new_class = c("mo_uncertainties", "data.frame"))
}

#' @rdname as.mo
#' @export
mo_renamed <- function() {
  add_MO_lookup_to_AMR_env()
  x <- AMR_env$mo_renamed

  x$new <- synonym_mo_to_accepted_mo(x$old)
  mo_old <- AMR_env$MO_lookup$fullname[match(x$old, AMR_env$MO_lookup$mo)]
  mo_new <- AMR_env$MO_lookup$fullname[match(x$new, AMR_env$MO_lookup$mo)]
  ref_old <- AMR_env$MO_lookup$ref[match(x$old, AMR_env$MO_lookup$mo)]
  ref_new <- AMR_env$MO_lookup$ref[match(x$new, AMR_env$MO_lookup$mo)]

  df_renamed <- data.frame(
    old = mo_old,
    new = mo_new,
    ref_old = ref_old,
    ref_new = ref_new,
    stringsAsFactors = FALSE
  )
  df_renamed <- unique(df_renamed)
  df_renamed <- df_renamed[order(df_renamed$old), , drop = FALSE]
  set_clean_class(df_renamed, new_class = c("mo_renamed", "data.frame"))
}

#' @rdname as.mo
#' @export
mo_failures <- function() {
  AMR_env$mo_failures
}

#' @rdname as.mo
#' @export
mo_reset_session <- function() {
  if (NROW(AMR_env$mo_previously_coerced) > 0) {
    message_("Reset ", nr2char(NROW(AMR_env$mo_previously_coerced)), " previously matched input value", ifelse(NROW(AMR_env$mo_previously_coerced) > 1, "s", ""), ".")
    AMR_env$mo_previously_coerced <- AMR_env$mo_previously_coerced[0, , drop = FALSE]
    AMR_env$mo_uncertainties <- AMR_env$mo_uncertainties[0, , drop = FALSE]
  } else {
    message_("No previously matched input values to reset.")
  }
}

#' @rdname as.mo
#' @export
mo_cleaning_regex <- function() {
  parts_to_remove <- c(
    "e?spp([^a-z]+|$)", "e?ssp([^a-z]+|$)", "e?ss([^a-z]+|$)", "e?sp([^a-z]+|$)", "e?subsp", "sube?species", "e?species",
    "biovar[a-z]*", "biotype", "serovar[a-z]*", "var([^a-z]+|$)", "serogr.?up[a-z]*",
    "titer", "dummy", "Ig[ADEGM]", " ?[a-z-]+[-](resistant|susceptible) ?"
  )

  paste0(
    "(",
    "[^A-Za-z- \\(\\)\\[\\]{}]+",
    "|",
    "([({]|\\[).+([})]|\\])",
    "|(^| )(",
    paste0(parts_to_remove[order(1 - nchar(parts_to_remove))], collapse = "|"),
    "))"
  )
}

# UNDOCUMENTED METHODS ----------------------------------------------------

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(pillar::pillar_shaft, mo)
pillar_shaft.mo <- function(x, ...) {
  add_MO_lookup_to_AMR_env()
  out <- trimws(format(x))
  # grey out the kingdom (part until first "_")
  out[!is.na(x)] <- gsub("^([A-Z]+_)(.*)", paste0(font_subtle("\\1"), "\\2"), out[!is.na(x)], perl = TRUE)
  # and grey out every _
  out[!is.na(x)] <- gsub("_", font_subtle("_"), out[!is.na(x)])

  # markup NA and UNKNOWN
  out[is.na(x)] <- font_na("  NA")
  out[x == "UNKNOWN"] <- font_na("  UNKNOWN")

  # markup manual codes
  out[x %in% AMR_env$MO_lookup$mo & !x %in% AMR::microorganisms$mo] <- font_blue(out[x %in% AMR_env$MO_lookup$mo & !x %in% AMR::microorganisms$mo], collapse = NULL)

  df <- tryCatch(get_current_data(arg_name = "x", call = 0),
    error = function(e) NULL
  )
  if (!is.null(df)) {
    mo_cols <- vapply(FUN.VALUE = logical(1), df, is.mo)
  } else {
    mo_cols <- NULL
  }

  all_mos <- c(AMR_env$MO_lookup$mo, NA)
  if (!all(x %in% all_mos) ||
    (!is.null(df) && !all(unlist(df[, which(mo_cols), drop = FALSE]) %in% all_mos))) {
    # markup old mo codes
    out[!x %in% all_mos] <- font_italic(
      font_na(x[!x %in% all_mos],
        collapse = NULL
      ),
      collapse = NULL
    )
    # throw a warning with the affected column name(s)
    if (!is.null(mo_cols)) {
      col <- paste0("Column ", vector_or(colnames(df)[mo_cols], quotes = TRUE, sort = FALSE))
    } else {
      col <- "The data"
    }
    warning_(
      col, " contains old MO codes (from a previous AMR package version). ",
      "Please update your MO codes with `as.mo()`.",
      call = FALSE
    )
  }

  # add the names to the bugs as mouse-over!
  if (tryCatch(isTRUE(getExportedValue("ansi_has_hyperlink_support", ns = asNamespace("cli"))()), error = function(e) FALSE)) {
    out[!x %in% c("UNKNOWN", NA)] <- font_url(
      url = paste0(
        x[!x %in% c("UNKNOWN", NA)], ": ",
        mo_name(x[!x %in% c("UNKNOWN", NA)], keep_synonyms = TRUE)
      ),
      txt = out[!x %in% c("UNKNOWN", NA)]
    )
  }

  # make it always fit exactly
  max_char <- max(nchar(x))
  if (is.na(max_char)) {
    max_char <- 12
  }
  create_pillar_column(out,
    align = "left",
    width = max_char + ifelse(any(x %in% c(NA, "UNKNOWN")), 2, 0)
  )
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(pillar::type_sum, mo)
type_sum.mo <- function(x, ...) {
  "mo"
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(cleaner::freq, mo)
freq.mo <- function(x, ...) {
  x_noNA <- as.mo(x[!is.na(x)]) # as.mo() to get the newest mo codes
  grams <- mo_gramstain(x_noNA, language = NULL)
  digits <- list(...)$digits
  if (is.null(digits)) {
    digits <- 2
  }
  cleaner::freq.default(
    x = x,
    ...,
    .add_header = list(
      `Gram-negative` = paste0(
        format(sum(grams == "Gram-negative", na.rm = TRUE),
          big.mark = " ",
          decimal.mark = "."
        ),
        " (", percentage(sum(grams == "Gram-negative", na.rm = TRUE) / length(grams),
          digits = digits
        ),
        ")"
      ),
      `Gram-positive` = paste0(
        format(sum(grams == "Gram-positive", na.rm = TRUE),
          big.mark = " ",
          decimal.mark = "."
        ),
        " (", percentage(sum(grams == "Gram-positive", na.rm = TRUE) / length(grams),
          digits = digits
        ),
        ")"
      ),
      `Nr. of genera` = pm_n_distinct(mo_genus(x_noNA, language = NULL)),
      `Nr. of species` = pm_n_distinct(paste(
        mo_genus(x_noNA, language = NULL),
        mo_species(x_noNA, language = NULL)
      ))
    )
  )
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(skimr::get_skimmers, mo)
get_skimmers.mo <- function(column) {
  mo <- as.mo(column, keep_synonyms = TRUE, language = NULL, info = FALSE)
  mo <- mo[!is.na(mo)]
  spp <- mo[mo_species(mo, keep_synonyms = TRUE, language = NULL, info = FALSE) != ""]
  skimr::sfl(
    skim_type = "mo",
    n_unique = ~ length(unique(mo)),
    gram_negative = ~ sum(mo_is_gram_negative(mo, keep_synonyms = TRUE, language = NULL, info = FALSE), na.rm = TRUE),
    gram_positive = ~ sum(mo_is_gram_positive(mo, keep_synonyms = TRUE, language = NULL, info = FALSE), na.rm = TRUE),
    yeast = ~ sum(mo_is_yeast(mo, keep_synonyms = TRUE, language = NULL, info = FALSE), na.rm = TRUE),
    top_genus = ~ names(sort(-table(mo_genus(mo, keep_synonyms = TRUE, language = NULL, info = FALSE))))[1L],
    top_species = ~ names(sort(-table(mo_name(spp, keep_synonyms = TRUE, language = NULL, info = FALSE))))[1L],
  )
}

#' @method print mo
#' @export
#' @noRd
print.mo <- function(x, print.shortnames = FALSE, ...) {
  add_MO_lookup_to_AMR_env()
  cat("Class 'mo'\n")
  x_names <- names(x)
  if (is.null(x_names) & print.shortnames == TRUE) {
    x_names <- tryCatch(mo_shortname(x, ...), error = function(e) NULL)
  }
  x <- as.character(x)
  names(x) <- x_names
  if (!all(x %in% c(AMR_env$MO_lookup$mo, NA))) {
    warning_(
      "Some MO codes are from a previous AMR package version. ",
      "Please update the MO codes with `as.mo()`.",
      call = FALSE
    )
  }
  print.default(x, quote = FALSE)
}

#' @method summary mo
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- object
  top_3 <- names(sort(-table(x[!is.na(x)])))[1:3]
  out <- c(
    "Class" = "mo",
    "<NA>" = length(x[is.na(x)]),
    "Unique" = length(unique(x[!is.na(x)])),
    "#1" = top_3[1],
    "#2" = top_3[2],
    "#3" = top_3[3]
  )
  class(out) <- c("summaryDefault", "table")
  out
}

#' @method as.data.frame mo
#' @export
#' @noRd
as.data.frame.mo <- function(x, ...) {
  add_MO_lookup_to_AMR_env()
  if (!all(x %in% c(AMR_env$MO_lookup$mo, NA))) {
    warning_(
      "The data contains old MO codes (from a previous AMR package version). ",
      "Please update your MO codes with `as.mo()`."
    )
  }
  nm <- deparse1(substitute(x))
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @method [ mo
#' @export
#' @noRd
"[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [[ mo
#' @export
#' @noRd
"[[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [<- mo
#' @export
#' @noRd
"[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  # must only contain valid MOs
  add_MO_lookup_to_AMR_env()
  return_after_integrity_check(y, "microorganism code", as.character(AMR_env$MO_lookup$mo))
}
#' @method [[<- mo
#' @export
#' @noRd
"[[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  # must only contain valid MOs
  add_MO_lookup_to_AMR_env()
  return_after_integrity_check(y, "microorganism code", as.character(AMR_env$MO_lookup$mo))
}
#' @method c mo
#' @export
#' @noRd
c.mo <- function(...) {
  x <- list(...)[[1L]]
  y <- NextMethod()
  attributes(y) <- attributes(x)
  add_MO_lookup_to_AMR_env()
  return_after_integrity_check(y, "microorganism code", as.character(AMR_env$MO_lookup$mo))
}

#' @method unique mo
#' @export
#' @noRd
unique.mo <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep mo
#' @export
#' @noRd
rep.mo <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method print mo_uncertainties
#' @export
#' @noRd
print.mo_uncertainties <- function(x, n = 10, ...) {
  more_than_50 <- FALSE
  if (NROW(x) == 0) {
    cat(word_wrap("No uncertainties to show. Only uncertainties of the last call to `as.mo()` or any `mo_*()` function are stored.\n\n", add_fn = font_blue))
    return(invisible(NULL))
  } else if (NROW(x) > 50) {
    more_than_50 <- TRUE
    x <- x[1:50, , drop = FALSE]
  }

  cat(word_wrap("Matching scores are based on the resemblance between the input and the full taxonomic name, and the pathogenicity in humans. See `?mo_matching_score`.\n\n", add_fn = font_blue))

  add_MO_lookup_to_AMR_env()

  col_red <- function(x) font_rose_bg(x, collapse = NULL)
  col_orange <- function(x) font_orange_bg(x, collapse = NULL)
  col_yellow <- function(x) font_yellow_bg(x, collapse = NULL)
  col_green <- function(x) font_green_bg(x, collapse = NULL)

  if (has_colour()) {
    cat(word_wrap("Colour keys: ",
      col_red(" 0.000-0.549 "),
      col_orange(" 0.550-0.649 "),
      col_yellow(" 0.650-0.749 "),
      col_green(" 0.750-1.000"),
      add_fn = font_blue
    ), font_green_bg(" "), "\n", sep = "")
  }

  score_set_colour <- function(text, scores) {
    # set colours to scores
    text[scores >= 0.75] <- col_green(text[scores >= 0.75])
    text[scores >= 0.65 & scores < 0.75] <- col_yellow(text[scores >= 0.65 & scores < 0.75])
    text[scores >= 0.55 & scores < 0.65] <- col_orange(text[scores >= 0.55 & scores < 0.65])
    text[scores < 0.55] <- col_red(text[scores < 0.55])
    text
  }

  txt <- ""
  any_maxed_out <- FALSE
  for (i in seq_len(nrow(x))) {
    if (x[i, ]$candidates != "") {
      candidates <- unlist(strsplit(x[i, ]$candidates, ", ", fixed = TRUE))
      if (length(candidates) > n) {
        any_maxed_out <- TRUE
        candidates <- candidates[seq_len(n)]
      }
      scores <- mo_matching_score(x = x[i, ]$input, n = candidates)
      n_candidates <- length(candidates)

      candidates_formatted <- italicise(candidates)
      scores_formatted <- trimws(formatC(round(scores, 3), format = "f", digits = 3))
      scores_formatted <- score_set_colour(scores_formatted, scores)

      # sort on descending scores
      candidates_formatted <- candidates_formatted[order(1 - scores)]
      scores_formatted <- scores_formatted[order(1 - scores)]

      candidates <- word_wrap(
        paste0(
          "Also matched: ",
          vector_and(
            paste0(
              candidates_formatted,
              font_blue(paste0(" (", scores_formatted, ")"), collapse = NULL)
            ),
            quotes = FALSE, sort = FALSE
          )
        ),
        extra_indent = nchar("Also matched: "),
        width = 0.9 * getOption("width", 100)
      )
    } else {
      candidates <- ""
    }

    score <- mo_matching_score(
      x = x[i, ]$input,
      n = x[i, ]$fullname
    )
    score_formatted <- trimws(formatC(round(score, 3), format = "f", digits = 3))
    txt <- paste(txt,
      paste0(
        paste0(
          "", strrep(font_grey("-"), times = getOption("width", 100)), "\n",
          '"', x[i, ]$original_input, '"',
          " -> ",
          paste0(
            font_bold(italicise(x[i, ]$fullname)),
            " (", x[i, ]$mo, ", ", score_set_colour(score_formatted, score), ")"
          )
        ),
        collapse = "\n"
      ),
      ifelse(x[i, ]$mo %in% AMR_env$MO_lookup$mo[which(AMR_env$MO_lookup$status == "synonym")],
        paste0(
          strrep(" ", nchar(x[i, ]$original_input) + 6),
          ifelse(x[i, ]$keep_synonyms == FALSE,
            # Add note if result was coerced to accepted taxonomic name
            font_red(paste0("This outdated taxonomic name was converted to ", font_italic(AMR_env$MO_lookup$fullname[match(synonym_mo_to_accepted_mo(x[i, ]$mo), AMR_env$MO_lookup$mo)], collapse = NULL), " (", synonym_mo_to_accepted_mo(x[i, ]$mo), ")."), collapse = NULL),
            # Or add note if result is currently another taxonomic name
            font_red(paste0(font_bold("Note: "), "The current name is ", font_italic(AMR_env$MO_lookup$fullname[match(synonym_mo_to_accepted_mo(x[i, ]$mo), AMR_env$MO_lookup$mo)], collapse = NULL), " (", AMR_env$MO_lookup$ref[match(synonym_mo_to_accepted_mo(x[i, ]$mo), AMR_env$MO_lookup$mo)], ")."), collapse = NULL)
          )
        ),
        ""
      ),
      candidates,
      sep = "\n"
    )
    txt <- gsub("[\n]+", "\n", txt)
    # remove first and last break
    txt <- gsub("(^[\n]|[\n]$)", "", txt)
    txt <- paste0("\n", txt, "\n")
  }

  cat(txt)
  if (isTRUE(any_maxed_out)) {
    cat(font_blue(word_wrap("\nOnly the first ", n, " other matches of each record are shown. Run `print(mo_uncertainties(), n = ...)` to view more entries, or save `mo_uncertainties()` to an object.")))
  }
  if (isTRUE(more_than_50)) {
    cat(font_blue(word_wrap("\nOnly the first 50 uncertainties are shown. Run `View(mo_uncertainties())` to view all entries, or save `mo_uncertainties()` to an object.")))
  }
}

#' @method print mo_renamed
#' @export
#' @noRd
print.mo_renamed <- function(x, extra_txt = "", n = 25, ...) {
  if (NROW(x) == 0) {
    cat(word_wrap("No renamed taxonomy to show. Only renamed taxonomy of the last call of `as.mo()` or any `mo_*()` function are stored.\n", add_fn = font_blue))
    return(invisible(NULL))
  }

  x$ref_old[!is.na(x$ref_old)] <- paste0(" (", gsub("et al.", font_italic("et al."), x$ref_old[!is.na(x$ref_old)], fixed = TRUE), ")")
  x$ref_new[!is.na(x$ref_new)] <- paste0(" (", gsub("et al.", font_italic("et al."), x$ref_new[!is.na(x$ref_new)], fixed = TRUE), ")")
  x$ref_old[is.na(x$ref_old)] <- " (author unknown)"
  x$ref_new[is.na(x$ref_new)] <- " (author unknown)"

  rows <- seq_len(min(NROW(x), n))

  message_(
    "The following microorganism", ifelse(NROW(x) > 1, "s were", " was"), " taxonomically renamed", extra_txt, ":\n",
    paste0("  ", AMR_env$bullet_icon, " ", font_italic(x$old[rows], collapse = NULL), x$ref_old[rows],
      "  ->  ", font_italic(x$new[rows], collapse = NULL), x$ref_new[rows],
      collapse = "\n"
    ),
    ifelse(NROW(x) > n, paste0("\n\nOnly the first ", n, " (out of ", NROW(x), ") are shown. Run `print(mo_renamed(), n = ...)` to view more entries (might be slow), or save `mo_renamed()` to an object."), "")
  )
}

# UNDOCUMENTED HELPER FUNCTIONS -------------------------------------------

convert_colloquial_input <- function(x) {
  x.bak <- trimws2(x)
  x <- trimws2(tolower(x))
  out <- rep(NA_character_, length(x))

  # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRPB)
  out[x %like_case% "^g[abcdefghijkl]s$"] <- gsub("g([abcdefghijkl])s",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "^g[abcdefghijkl]s$"],
    perl = TRUE
  )
  # Streptococci in different languages, like "estreptococos grupo B"
  out[x %like_case% "strepto[ck]o[ck][a-zA-Z ]* [abcdefghijkl]$"] <- gsub(".*e?strepto[ck]o[ck].* ([abcdefghijkl])$",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "strepto[ck]o[ck][a-zA-Z ]* [abcdefghijkl]$"],
    perl = TRUE
  )
  out[x %like_case% "strep[a-z]* group [abcdefghijkl]$"] <- gsub(".* ([abcdefghijkl])$",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "strep[a-z]* group [abcdefghijkl]$"],
    perl = TRUE
  )
  out[x %like_case% "group [abcdefghijkl] strepto[ck]o[ck]"] <- gsub(".*group ([abcdefghijkl]) strepto[ck]o[ck].*",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "group [abcdefghijkl] strepto[ck]o[ck]"],
    perl = TRUE
  )
  out[x %like_case% "ha?emoly.*strep"] <- "B_STRPT_HAEM"
  out[x %like_case% "(strepto.* [abcg, ]{2,4}$)"] <- "B_STRPT_ABCG"
  out[x %like_case% "(strepto.* mil+er+i|^mgs[^a-z]*$)"] <- "B_STRPT_MILL"
  out[x %like_case% "mil+er+i gr"] <- "B_STRPT_MILL"
  out[x %like_case% "((strepto|^s).* viridans|^vgs[^a-z]*$)"] <- "B_STRPT_VIRI"
  out[x %like_case% "(viridans.* (strepto|^s).*|^vgs[^a-z]*$)"] <- "B_STRPT_VIRI"
  out[x %like_case% "meningo[ck]o[ck](ken)?$"] <- "B_NESSR_MNNG"
  out[x %like_case% "pneumo[ck]o[ck](ken)?$"] <- "B_STRPT_PNMN"

  # Salmonella in different languages, like "Salmonella grupo B"
  out[x %like_case% "salmonella.* [abcdefgh]$"] <- gsub(".*salmonella.* ([abcdefgh])$",
    "B_SLMNL_GRP\\U\\1",
    x[x %like_case% "salmonella.* [abcdefgh]$"],
    perl = TRUE
  )
  out[x %like_case% "group [abcdefgh] salmonella"] <- gsub(".*group ([abcdefgh]) salmonella*",
    "B_SLMNL_GRP\\U\\1",
    x[x %like_case% "group [abcdefgh] salmonella"],
    perl = TRUE
  )

  # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese)
  out[x %like_case% "([ck]oagulas[ea].negatie?[vf]|^[ck]o?ns[^a-z]*$)"] <- "B_STPHY_CONS"
  out[x %like_case% "([ck]oagulas[ea].positie?[vf]|^[ck]o?ps[^a-z]*$)"] <- "B_STPHY_COPS"

  # Gram stains
  out[x %like_case% "gram[ -]?neg.*"] <- "B_GRAMN"
  out[x %like_case% "( |^)gram[-]( |$)"] <- "B_GRAMN"
  out[x %like_case% "gram[ -]?pos.*"] <- "B_GRAMP"
  out[x %like_case% "( |^)gram[+]( |$)"] <- "B_GRAMP"
  out[x %like_case% "anaerob[a-z]+ .*gram[ -]?neg.*"] <- "B_ANAER-NEG"
  out[x %like_case% "anaerob[a-z]+ .*gram[ -]?pos.*"] <- "B_ANAER-POS"
  out[is.na(out) & x %like_case% "anaerob[a-z]+ (micro)?.*organism"] <- "B_ANAER"
  out[is.na(out) & x %like_case% "anaerob[a-z]+ bacter"] <- "B_ANAER"

  # coryneform bacteria
  out[x %like_case% "^coryneform"] <- "B_CORYNF"

  # yeasts and fungi
  out[x %like_case% "(^| )yeast?"] <- "F_YEAST"
  out[x %like_case% "(^| )fung(us|i)"] <- "F_FUNGUS"

  # protozoa
  out[x %like_case% "protozo"] <- "P_PROTOZOAN" # to hit it with most languages, and "protozo" does not occur in the microorganisms data set for anything else

  # trivial names known to the field
  out[x %like_case% "meningo[ck]o[ck]"] <- "B_NESSR_MNNG"
  out[x %like_case% "gono[ck]o[ck]"] <- "B_NESSR_GNRR"
  out[x %like_case% "pneumo[ck]o[ck]"] <- "B_STRPT_PNMN"
  out[x %like_case% "hacek"] <- "B_HACEK"
  out[x %like_case% "haemophilus" & x %like_case% "aggregatibacter" & x %like_case% "cardiobacterium" & x %like_case% "eikenella" & x %like_case% "kingella"] <- "B_HACEK"
  out[x %like_case% "slow.* grow.* mycobact"] <- "B_MYCBC_SGM"
  out[x %like_case% "rapid.* grow.* mycobact"] <- "B_MYCBC_RGM"

  # unexisting names (con is the WHONET code for contamination)
  out[x %in% c("con", "other", "none", "unknown") | x %like_case% "virus"] <- "UNKNOWN"

  # WHONET has a lot of E. coli and Vibrio cholerae names
  out[x %like_case% "escherichia coli"] <- "B_ESCHR_COLI"
  out[x %like_case% "vibrio cholerae"] <- "B_VIBRI_CHLR"

  out
}

italicise <- function(x) {
  if (!has_colour()) {
    return(x)
  }
  out <- font_italic(x, collapse = NULL)
  # city-like serovars of Salmonella (start with a capital)
  out[x %like_case% "Salmonella [A-Z]"] <- paste(
    font_italic("Salmonella"),
    gsub("Salmonella ", "", x[x %like_case% "Salmonella [A-Z]"])
  )
  # streptococcal groups
  out[x %like_case% "Streptococcus [A-Z]"] <- paste(
    font_italic("Streptococcus"),
    gsub("Streptococcus ", "", x[x %like_case% "Streptococcus [A-Z]"])
  )
  # be sure not to make these italic
  out <- gsub("([ -]*)(Group|group|Complex|complex)(\033\\[23m)?", "\033[23m\\1\\2", out, perl = TRUE)
  out <- gsub("(\033\\[3m)?(Beta[-]haemolytic|Coagulase[-](postive|negative)) ", "\\2 \033[3m", out, perl = TRUE)
  out
}

nr2char <- function(x) {
  if (x %in% c(1:10)) {
    v <- c(
      "one" = 1, "two" = 2, "three" = 3, "four" = 4, "five" = 5,
      "six" = 6, "seven" = 7, "eight" = 8, "nine" = 9, "ten" = 10
    )
    names(v[x])
  } else {
    x
  }
}

parse_and_convert <- function(x) {
  if (tryCatch(is.character(x) && all(Encoding(x) == "unknown", na.rm = TRUE), error = function(e) FALSE)) {
    out <- x
  } else {
    out <- tryCatch(
      {
        if (!is.null(dim(x))) {
          if (NCOL(x) > 2) {
            stop("a maximum of two columns is allowed", call. = FALSE)
          } else if (NCOL(x) == 2) {
            # support Tidyverse selection like: df %>% select(colA, colB)
            # paste these columns together
            x <- as.data.frame(x, stringsAsFactors = FALSE)
            colnames(x) <- c("A", "B")
            x <- paste(x$A, x$B)
          } else {
            # support Tidyverse selection like: df %>% select(colA)
            x <- as.data.frame(x, stringsAsFactors = FALSE)[[1]]
          }
        }
        parsed <- iconv(as.character(x), to = "UTF-8")
        parsed[is.na(parsed) & !is.na(x)] <- iconv(x[is.na(parsed) & !is.na(x)], from = "Latin1", to = "ASCII//TRANSLIT")
        parsed <- gsub('"', "", parsed, fixed = TRUE)
        parsed
      },
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    ) # this will also be thrown when running `as.mo(no_existing_object)`
  }
  out <- trimws2(out)
  out <- gsub(" +", " ", out, perl = TRUE)
  out <- gsub(" ?/ ? ", "/", out, perl = TRUE)
  out
}

replace_old_mo_codes <- function(x, property) {
  # this function transform old MO codes to current codes, such as:
  # B_ESCH_COL (AMR v0.5.0) -> B_ESCHR_COLI
  ind <- x %like_case% "^[A-Z]_[A-Z_]+$" & !x %in% AMR_env$MO_lookup$mo
  if (any(ind, na.rm = TRUE)) {
    add_MO_lookup_to_AMR_env()
    # get the ones that match
    affected <- x[ind]
    affected_unique <- unique(affected)
    all_direct_matches <- TRUE
    # find their new codes, once per code
    solved_unique <- unlist(lapply(
      strsplit(affected_unique, ""),
      function(m) {
        kingdom <- paste0("^", m[1])
        name <- m[3:length(m)]
        name[name == "_"] <- " "
        name <- tolower(paste0(name, ".*", collapse = ""))
        name <- gsub(" .*", " ", name, fixed = TRUE)
        name <- paste0("^", name)
        results <- AMR_env$MO_lookup$mo[AMR_env$MO_lookup$kingdom %like_case% kingdom &
          AMR_env$MO_lookup$fullname_lower %like_case% name]
        if (length(results) > 1) {
          all_direct_matches <<- FALSE
        }
        results[1L]
      }
    ), use.names = FALSE)
    solved <- solved_unique[match(affected, affected_unique)]
    # assign on places where a match was found
    x[ind] <- solved
    n_matched <- length(affected[!is.na(affected)])
    n_solved <- length(affected[!is.na(solved)])
    n_unsolved <- length(affected[is.na(solved)])
    n_unique <- length(affected_unique[!is.na(affected_unique)])
    if (n_unique < n_matched) {
      n_unique <- paste0(n_unique, " unique, ")
    } else {
      n_unique <- ""
    }
    if (property != "mo") {
      warning_(
        "in `mo_", property, "()`: the input contained ", n_matched,
        " old MO code", ifelse(n_matched == 1, "", "s"),
        " (", n_unique, "from a previous AMR package version). ",
        "Please update your MO codes with `as.mo()` to increase speed."
      )
    } else {
      warning_(
        "in `as.mo()`: the input contained ", n_matched,
        " old MO code", ifelse(n_matched == 1, "", "s"),
        " (", n_unique, "from a previous AMR package version). ",
        n_solved, " old MO code", ifelse(n_solved == 1, "", "s"),
        ifelse(n_solved == 1, " was", " were"),
        ifelse(all_direct_matches, " updated ", font_bold(" guessed ")),
        "to ", ifelse(n_solved == 1, "a ", ""),
        "currently used MO code", ifelse(n_solved == 1, "", "s"),
        ifelse(n_unsolved > 0,
          paste0(" and ", n_unsolved, " old MO code", ifelse(n_unsolved == 1, "", "s"), " could not be updated."),
          "."
        )
      )
    }
  }
  x
}

replace_ignore_pattern <- function(x, ignore_pattern) {
  if (!is.null(ignore_pattern) && !identical(trimws2(ignore_pattern), "")) {
    ignore_cases <- x %like% ignore_pattern
    if (sum(ignore_cases) > 0) {
      message_(
        "The following input was ignored by `ignore_pattern = \"", ignore_pattern, "\"`: ",
        vector_and(x[ignore_cases], quotes = TRUE)
      )
      x[ignore_cases] <- NA_character_
    }
  }
  x
}

repair_reference_df <- function(reference_df) {
  if (is.null(reference_df)) {
    return(NULL)
  }
  # has valid own reference_df
  reference_df <- reference_df %pm>%
    pm_filter(!is.na(mo))

  # keep only first two columns, second must be mo
  if (colnames(reference_df)[1] == "mo") {
    reference_df <- reference_df %pm>% pm_select(2, "mo")
  } else {
    reference_df <- reference_df %pm>% pm_select(1, "mo")
  }

  # remove factors, just keep characters
  colnames(reference_df)[1] <- "x"
  reference_df[, "x"] <- as.character(reference_df[, "x", drop = TRUE])
  reference_df[, "mo"] <- as.character(reference_df[, "mo", drop = TRUE])

  # some MO codes might be old
  reference_df[, "mo"] <- as.mo(reference_df[, "mo", drop = TRUE], reference_df = NULL)
  reference_df
}

get_mo_uncertainties <- function() {
  remember <- list(
    uncertainties = AMR_env$mo_uncertainties,
    failures = AMR_env$mo_failures
  )
  # empty them, otherwise e.g. mo_shortname("Chlamydophila psittaci") will give 3 notes
  AMR_env$mo_uncertainties <- NULL
  AMR_env$mo_failures <- NULL
  remember
}

load_mo_uncertainties <- function(metadata) {
  AMR_env$mo_uncertainties <- metadata$uncertainties
  AMR_env$mo_failures <- metadata$failures
}

synonym_mo_to_accepted_mo <- function(x, fill_in_accepted = FALSE, dataset = AMR_env$MO_lookup) {
  # `dataset` is an argument so that it can be used in the regeneration of the microorganisms data set
  if (identical(dataset, AMR_env$MO_lookup)) {
    add_MO_lookup_to_AMR_env()
    dataset <- AMR_env$MO_lookup
  }

  out <- x
  is_still_synonym <- dataset$status[match(out, dataset$mo)] == "synonym"
  limit <- 0
  while (any(is_still_synonym, na.rm = TRUE) && limit < 5) {
    limit <- limit + 1

    # make sure to get the latest name, e.g. Fusarium pulicaris robiniae was first renamed to Fusarium roseum, then to Fusarium sambucinum
    # we need the MO of Fusarium pulicaris robiniae to return the MO of Fusarium sambucinum
    must_be_corrected <- !is.na(is_still_synonym) & is_still_synonym
    x_gbif <- dataset$gbif_renamed_to[match(out, dataset$mo)]
    x_mycobank <- dataset$mycobank_renamed_to[match(out, dataset$mo)]
    x_lpsn <- dataset$lpsn_renamed_to[match(out, dataset$mo)]

    out[must_be_corrected & !is.na(x_gbif)] <- dataset$mo[match(x_gbif[must_be_corrected & !is.na(x_gbif)], dataset$gbif)]
    out[must_be_corrected & !is.na(x_mycobank)] <- dataset$mo[match(x_mycobank[must_be_corrected & !is.na(x_mycobank)], dataset$mycobank)]
    out[must_be_corrected & !is.na(x_lpsn)] <- dataset$mo[match(x_lpsn[must_be_corrected & !is.na(x_lpsn)], dataset$lpsn)]

    is_still_synonym <- dataset$status[match(out, dataset$mo)] == "synonym"
  }

  x_no_synonym <- dataset$status[match(x, dataset$mo)] != "synonym"
  out[x_no_synonym] <- NA_character_
  if (isTRUE(fill_in_accepted)) {
    out[!is.na(x_no_synonym) & x_no_synonym] <- x[!is.na(x_no_synonym) & x_no_synonym]
  }

  out[is.na(match(x, dataset$mo))] <- NA_character_
  out
}

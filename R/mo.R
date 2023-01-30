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

#' Transform Input to a Microorganism Code
#'
#' Use this function to determine a valid microorganism code ([`mo`]). Determination is done using intelligent rules and the complete taxonomic kingdoms `r vector_and(unique(microorganisms$kingdom[which(!grepl("(unknown|Fungi)", microorganisms$kingdom))]), quotes = FALSE)`, and most microbial species from the kingdom Fungi (see *Source*). The input can be almost anything: a full name (like `"Staphylococcus aureus"`), an abbreviated name (such as `"S. aureus"`), an abbreviation known in the field (such as `"MRSA"`), or just a genus. See *Examples*.
#' @param x a [character] vector or a [data.frame] with one or two columns
#' @param Becker a [logical] to indicate whether staphylococci should be categorised into coagulase-negative staphylococci ("CoNS") and coagulase-positive staphylococci ("CoPS") instead of their own species, according to Karsten Becker *et al.* (see Source).
#'
#' This excludes *Staphylococcus aureus* at default, use `Becker = "all"` to also categorise *S. aureus* as "CoPS".
#' @param Lancefield a [logical] to indicate whether a beta-haemolytic *Streptococcus* should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield (see Source). These streptococci will be categorised in their first group, e.g. *Streptococcus dysgalactiae* will be group C, although officially it was also categorised into groups G and L.
#'
#' This excludes enterococci at default (who are in group D), use `Lancefield = "all"` to also categorise all enterococci as group D.
#' @param minimum_matching_score a numeric value to set as the lower limit for the [MO matching score][mo_matching_score()]. When left blank, this will be determined automatically based on the character length of `x`, its [taxonomic kingdom][microorganisms] and [human pathogenicity][mo_matching_score()].
#' @param keep_synonyms a [logical] to indicate if old, previously valid taxonomic names must be preserved and not be corrected to currently accepted names. The default is `FALSE`, which will return a note if old taxonomic names were processed. The default can be set with the option [`AMR_keep_synonyms`][AMR-options], i.e. `options(AMR_keep_synonyms = TRUE)` or `options(AMR_keep_synonyms = FALSE)`.
#' @param reference_df a [data.frame] to be used for extra reference when translating `x` to a valid [`mo`]. See [set_mo_source()] and [get_mo_source()] to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ignore_pattern a [regular expression][base::regex] (case-insensitive) of which all matches in `x` must return `NA`. This can be convenient to exclude known non-relevant input and can also be set with the option [`AMR_ignore_pattern`][AMR-options], e.g. `options(AMR_ignore_pattern = "(not reported|contaminated flora)")`.
#' @param remove_from_input a [regular expression][base::regex] (case-insensitive) to clean the input of `x`. Everything matched in `x` will be removed. At default, this is the outcome of [mo_cleaning_regex()], which removes texts between brackets and texts such as "species" and "serovar".
#' @param language language to translate text like "no growth", which defaults to the system language (see [get_AMR_locale()])
#' @param info a [logical] to indicate if a progress bar should be printed if more than 25 items are to be coerced, defaults to `TRUE` only in interactive mode
#' @param ... other arguments passed on to functions
#' @rdname as.mo
#' @aliases mo
#' @details
#' A microorganism (MO) code from this package (class: [`mo`]) is human readable and typically looks like these examples:
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
#'   \----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                             F (Fungi), PL (Plantae), P (Protozoa)
#' ```
#'
#' Values that cannot be coerced will be considered 'unknown' and will be returned as the MO code `UNKNOWN` with a warning.
#'
#' Use the [`mo_*`][mo_property()] functions to get properties based on the returned code, see *Examples*.
#'
#' The [as.mo()] function uses a novel [matching score algorithm][mo_matching_score()] (see *Matching Score for Microorganisms* below) to match input against the [available microbial taxonomy][microorganisms] in this package. This will lead to the effect that e.g. `"E. coli"` (a microorganism highly prevalent in humans) will return the microbial ID of *Escherichia coli* and not *Entamoeba coli* (a microorganism less prevalent in humans), although the latter would alphabetically come first. The algorithm uses data from the List of Prokaryotic names with Standing in Nomenclature (LPSN) and the Global Biodiversity Information Facility (GBIF) (see [microorganisms]).
#'
#' ### Coping with Uncertain Results
#'
#' Results of non-exact taxonomic input are based on their [matching score][mo_matching_score()]. The lowest allowed score can be set with the `minimum_matching_score` argument. At default this will be determined based on the character length of the input, and the [taxonomic kingdom][microorganisms] and [human pathogenicity][mo_matching_score()] of the taxonomic outcome. If values are matched with uncertainty, a message will be shown to suggest the user to evaluate the results with [mo_uncertainties()], which returns a [data.frame] with all specifications.
#'
#' To increase the quality of matching, the `remove_from_input` argument can be used to clean the input (i.e., `x`). This must be a [regular expression][base::regex] that matches parts of the input that should be removed before the input is matched against the [available microbial taxonomy][microorganisms]. It will be matched Perl-compatible and case-insensitive. The default value of `remove_from_input` is the outcome of the helper function [mo_cleaning_regex()].
#'
#' There are three helper functions that can be run after using the [as.mo()] function:
#' - Use [mo_uncertainties()] to get a [data.frame] that prints in a pretty format with all taxonomic names that were guessed. The output contains the matching score for all matches (see *Matching Score for Microorganisms* below).
#' - Use [mo_failures()] to get a [character] [vector] with all values that could not be coerced to a valid value.
#' - Use [mo_renamed()] to get a [data.frame] with all values that could be coerced based on old, previously accepted taxonomic names.
#'
#' ### Microbial Prevalence of Pathogens in Humans
#'
#' The coercion rules consider the prevalence of microorganisms in humans, which is available as the `prevalence` column in the [microorganisms] data set. The grouping into human pathogenic prevalence is explained in the section *Matching Score for Microorganisms* below.
#' @inheritSection mo_matching_score Matching Score for Microorganisms
#'
#  (source as a section here, so it can be inherited by other man pages)
#' @section Source:
#' 1. Berends MS *et al.* (2022). **AMR: An R Package for Working with Antimicrobial Resistance Data**. *Journal of Statistical Software*, 104(3), 1-31; \doi{10.18637/jss.v104.i03}
#' 2. Becker K *et al.* (2014). **Coagulase-Negative Staphylococci.** *Clin Microbiol Rev.* 27(4): 870-926; \doi{10.1128/CMR.00109-13}
#' 3. Becker K *et al.* (2019). **Implications of identifying the recently defined members of the *S. aureus* complex, *S. argenteus* and *S. schweitzeri*: A position paper of members of the ESCMID Study Group for staphylococci and Staphylococcal Diseases (ESGS).** *Clin Microbiol Infect*; \doi{10.1016/j.cmi.2019.02.028}
#' 4. Becker K *et al.* (2020). **Emergence of coagulase-negative staphylococci** *Expert Rev Anti Infect Ther.* 18(4):349-366; \doi{10.1080/14787210.2020.1730813}
#' 5. Lancefield RC (1933). **A serological differentiation of human and other groups of hemolytic streptococci**. *J Exp Med.* 57(4): 571-95; \doi{10.1084/jem.57.4.571}
#' 6. Berends MS *et al.* (2022). **Trends in Occurrence and Phenotypic Resistance of Coagulase-Negative Staphylococci (CoNS) Found in Human Blood in the Northern Netherlands between 2013 and 2019** *Microorganisms* 10(9), 1801; \doi{10.3390/microorganisms10091801}
#' 7. `r TAXONOMY_VERSION$LPSN$citation` Accessed from <`r TAXONOMY_VERSION$LPSN$url`> on `r documentation_date(TAXONOMY_VERSION$LPSN$accessed_date)`.
#' 8. `r TAXONOMY_VERSION$GBIF$citation` Accessed from <`r TAXONOMY_VERSION$GBIF$url`> on `r documentation_date(TAXONOMY_VERSION$GBIF$accessed_date)`.
#' 9. `r TAXONOMY_VERSION$SNOMED$citation` URL: <`r TAXONOMY_VERSION$SNOMED$url`>
#' 10. Bartlett A *et al.* (2022). **A comprehensive list of bacterial pathogens infecting humans** *Microbiology* 168:001269; \doi{10.1099/mic.0.001269}
#' @export
#' @return A [character] [vector] with additional class [`mo`]
#' @seealso [microorganisms] for the [data.frame] that is being used to determine ID's.
#'
#' The [`mo_*`][mo_property()] functions (such as [mo_genus()], [mo_gramstain()]) to get properties based on the returned code.
#' @inheritSection AMR Reference Data Publicly Available
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
                  remove_from_input = mo_cleaning_regex(),
                  language = get_AMR_locale(),
                  info = interactive(),
                  ...) {
  meet_criteria(x, allow_class = c("mo", "data.frame", "list", "character", "numeric", "integer", "factor"), allow_NA = TRUE)
  meet_criteria(Becker, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(Lancefield, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)
  meet_criteria(minimum_matching_score, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE)
  meet_criteria(reference_df, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(ignore_pattern, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  language <- validate_language(language)
  meet_criteria(info, allow_class = "logical", has_length = 1)

  add_MO_lookup_to_AMR_env()

  if (tryCatch(all(x %in% c(AMR_env$MO_lookup$mo, NA)) &&
    isFALSE(Becker) &&
    isFALSE(Lancefield), error = function(e) FALSE)) {
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

  complexes <- x[trimws2(x_lower) %like_case% " (complex|group)$"]
  if (length(complexes) > 0 && identical(remove_from_input, mo_cleaning_regex()) && !any(AMR_env$MO_lookup$fullname[which(AMR_env$MO_lookup$source == "Added by user")] %like% "(group|complex)", na.rm = TRUE)) {
    warning_("in `as.mo()`: 'complex' and 'group' were ignored from the input in ", length(complexes), " case", ifelse(length(complexes) > 1, "s", ""), ", as they are currently not supported.\nYou can add your own microorganism with `add_custom_microorganisms()`.", call = FALSE)
  }

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
  out[is.na(out) & toupper(x) %in% AMR::microorganisms.codes$code] <- AMR::microorganisms.codes$mo[match(toupper(x)[is.na(out) & toupper(x) %in% AMR::microorganisms.codes$code], AMR::microorganisms.codes$code)]
  # From SNOMED ----
  if (any(is.na(out) & !is.na(x)) && any(is.na(out) & x %in% unlist(AMR_env$MO_lookup$snomed), na.rm = TRUE)) {
    # found this extremely fast gem here: https://stackoverflow.com/a/11002456/4575331
    out[is.na(out) & x %in% unlist(AMR_env$MO_lookup$snomed)] <- AMR_env$MO_lookup$mo[rep(seq_along(AMR_env$MO_lookup$snomed), vapply(FUN.VALUE = double(1), AMR_env$MO_lookup$snomed, length))[match(x[is.na(out) & x %in% unlist(AMR_env$MO_lookup$snomed)], unlist(AMR_env$MO_lookup$snomed))]]
  }
  # From other familiar output ----
  # such as Salmonella groups, colloquial names, etc.
  out[is.na(out)] <- convert_colloquial_input(x[is.na(out)])
  # From previous hits in this session ----
  old <- out
  out[is.na(out) & paste(x, minimum_matching_score) %in% AMR_env$mo_previously_coerced$x] <- AMR_env$mo_previously_coerced$mo[match(paste(x, minimum_matching_score)[is.na(out) & paste(x, minimum_matching_score) %in% AMR_env$mo_previously_coerced$x], AMR_env$mo_previously_coerced$x)]
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
    x[trimws2(x) %like% translate_into_language("no .*growth", language = language)] <- NA_character_
    x[trimws2(x) %like% paste0("^(", translate_into_language("no|not", language = language), ") ")] <- NA_character_

    # groups are in our taxonomic table with a capital G
    x <- gsub(" group( |$)", " Group\\1", x, perl = TRUE)

    # run over all unique leftovers
    x_unique <- unique(x[is.na(out) & !is.na(x)])

    # set up progress bar
    progress <- progress_ticker(n = length(x_unique), n_min = 10, print = info)
    on.exit(close(progress))

    msg <- character(0)

    # run it
    x_coerced <- vapply(FUN.VALUE = character(1), x_unique, function(x_search) {
      progress$tick()

      # some required cleaning steps
      x_out <- trimws2(x_search)
      # this applies the `remove_from_input` argument, which defaults to mo_cleaning_regex()
      x_out <- gsub(remove_from_input, " ", x_out, ignore.case = TRUE, perl = TRUE)
      x_out <- trimws2(gsub(" +", " ", x_out, perl = TRUE))
      x_search_cleaned <- x_out
      x_out <- tolower(x_out)

      # first check if cleaning led to an exact result, case-insensitive
      if (x_out %in% AMR_env$MO_lookup$fullname_lower) {
        return(as.character(AMR_env$MO_lookup$mo[match(x_out, AMR_env$MO_lookup$fullname_lower)]))
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
        filtr <- which(AMR_env$MO_lookup$full_first == substr(x_parts[1], 1, 1) & (AMR_env$MO_lookup$species_first == substr(x_parts[2], 1, 1) | AMR_env$MO_lookup$subspecies_first == substr(x_parts[2], 1, 1)))
      } else if (length(x_parts) > 3) {
        first_chars <- paste0("(^| )", "[", paste(substr(x_parts, 1, 1), collapse = ""), "]")
        filtr <- which(AMR_env$MO_lookup$full_first %like_case% first_chars)
      } else if (nchar(x_out) == 4) {
        # no space and 4 characters - probably a code such as STAU or ESCO
        msg <- c(msg, paste0("Input \"", x_search, "\" was assumed to be a microorganism code - tried to match on ", vector_and(c(substr(x_out, 1, 2), substr(x_out, 3, 4)), sort = FALSE)))
        filtr <- which(AMR_env$MO_lookup$fullname_lower %like_case% paste0("(^| )", substr(x_out, 1, 2), ".* ", substr(x_out, 3, 4)))
      } else if (nchar(x_out) <= 6) {
        # no space and 5-6 characters - probably a code such as STAAUR or ESCCOL
        first_part <- paste0(substr(x_out, 1, 2), "[a-z]*", substr(x_out, 3, 3))
        second_part <- substr(x_out, 4, nchar(x_out))
        msg <- c(msg, paste0("Input \"", x_search, "\" was assumed to be a microorganism code - tried to match on ", vector_and(c(gsub("[a-z]*", "(...)", first_part, fixed = TRUE), second_part), sort = FALSE)))
        filtr <- which(AMR_env$MO_lookup$fullname_lower %like_case% paste0("(^| )", first_part, ".* ", second_part))
      } else {
        filtr <- which(AMR_env$MO_lookup$full_first == substr(x_out, 1, 1))
      }
      if (length(filtr) == 0) {
        mo_to_search <- AMR_env$MO_lookup$fullname
      } else {
        mo_to_search <- AMR_env$MO_lookup$fullname[filtr]
      }

      AMR_env$mo_to_search <- mo_to_search
      # determine the matching score on the original search value
      m <- mo_matching_score(x = x_search_cleaned, n = mo_to_search)
      if (is.null(minimum_matching_score)) {
        minimum_matching_score_current <- min(0.6, min(10, nchar(x_search_cleaned)) * 0.08)
        # correct back for prevalence
        minimum_matching_score_current <- minimum_matching_score_current / AMR_env$MO_lookup$prevalence[match(mo_to_search, AMR_env$MO_lookup$fullname)]
        # correct back for kingdom
        minimum_matching_score_current <- minimum_matching_score_current / AMR_env$MO_lookup$kingdom_index[match(mo_to_search, AMR_env$MO_lookup$fullname)]
      } else {
        minimum_matching_score_current <- minimum_matching_score
      }

      if (sum(m >= minimum_matching_score_current) > 10) {
        # at least 10 are left over, make the ones under `m` NA
        m[m < minimum_matching_score_current] <- NA_real_
      }

      top_hits <- mo_to_search[order(m, decreasing = TRUE, na.last = NA)] # na.last = NA will remove the NAs
      if (length(top_hits) == 0) {
        warning_("No hits found for \"", x_search, "\" with minimum_matching_score = ", ifelse(is.null(minimum_matching_score), paste0("NULL (=", round(min(minimum_matching_score_current, na.rm = TRUE), 3), ")"), minimum_matching_score), ". Try setting this value lower or even to 0.", call = FALSE)
        result_mo <- NA_character_
      } else {
        result_mo <- AMR_env$MO_lookup$mo[match(top_hits[1], AMR_env$MO_lookup$fullname)]
        AMR_env$mo_uncertainties <- rbind(AMR_env$mo_uncertainties,
          data.frame(
            original_input = x_search,
            input = x_search_cleaned,
            fullname = top_hits[1],
            mo = result_mo,
            candidates = ifelse(length(top_hits) > 1, paste(top_hits[2:min(26, length(top_hits))], collapse = ", "), ""),
            minimum_matching_score = ifelse(is.null(minimum_matching_score), "NULL", minimum_matching_score),
            keep_synonyms = keep_synonyms,
            stringsAsFactors = FALSE
          ),
          stringsAsFactors = FALSE
        )
        # save to package env to save time for next time
        AMR_env$mo_previously_coerced <- unique(rbind(AMR_env$mo_previously_coerced,
          data.frame(
            x = paste(x_search, minimum_matching_score),
            mo = result_mo,
            stringsAsFactors = FALSE
          ),
          stringsAsFactors = FALSE
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
  lpsn_matches <- AMR_env$MO_lookup$lpsn_renamed_to[match(out, AMR_env$MO_lookup$mo)]
  lpsn_matches[!lpsn_matches %in% AMR_env$MO_lookup$lpsn] <- NA
  # GBIF only for non-bacteria, since we use LPSN as primary source for bacteria
  # (an example is Strep anginosus, renamed according to GBIF, not according to LPSN)
  gbif_matches <- AMR_env$MO_lookup$gbif_renamed_to[AMR_env$MO_lookup$kingdom != "Bacteria"][match(out, AMR_env$MO_lookup$mo[AMR_env$MO_lookup$kingdom != "Bacteria"])]
  gbif_matches[!gbif_matches %in% AMR_env$MO_lookup$gbif] <- NA
  AMR_env$mo_renamed <- list(
    old = out[!is.na(gbif_matches) | !is.na(lpsn_matches)],
    gbif_matches = gbif_matches[!is.na(gbif_matches) | !is.na(lpsn_matches)],
    lpsn_matches = lpsn_matches[!is.na(gbif_matches) | !is.na(lpsn_matches)]
  )
  if (isFALSE(keep_synonyms)) {
    out[which(!is.na(gbif_matches))] <- AMR_env$MO_lookup$mo[match(gbif_matches[which(!is.na(gbif_matches))], AMR_env$MO_lookup$gbif)]
    out[which(!is.na(lpsn_matches))] <- AMR_env$MO_lookup$mo[match(lpsn_matches[which(!is.na(lpsn_matches))], AMR_env$MO_lookup$lpsn)]
    if (isTRUE(info) && length(AMR_env$mo_renamed$old) > 0) {
      print(mo_renamed(), extra_txt = " (use `keep_synonyms = TRUE` to leave uncorrected)")
    }
  } else if (is.null(getOption("AMR_keep_synonyms")) && length(AMR_env$mo_renamed$old) > 0 && message_not_thrown_before("as.mo", "keep_synonyms_warning", entire_session = TRUE)) {
    # keep synonyms is TRUE, so check if any do have synonyms
    warning_("Function `as.mo()` returned ", nr2char(length(unique(AMR_env$mo_renamed$old))), " old taxonomic name", ifelse(length(unique(AMR_env$mo_renamed$old)) > 1, "s", ""), ". Use `as.mo(..., keep_synonyms = FALSE)` to clean the input to currently accepted taxonomic names, or set the R option `AMR_keep_synonyms` to `FALSE`. This warning will be shown once per session.", call = FALSE)
  }

  # Apply Becker ----
  if (isTRUE(Becker) || Becker == "all") {
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

    # 'MO_CONS' and 'MO_COPS' are 'mo' vectors created in R/_pre_commit_hook.R
    out[out %in% MO_CONS] <- "B_STPHY_CONS"
    out[out %in% MO_COPS] <- "B_STPHY_COPS"
    if (Becker == "all") {
      out[out == "B_STPHY_AURS"] <- "B_STPHY_COPS"
    }
  }

  # Apply Lancefield ----
  if (isTRUE(Lancefield) || Lancefield == "all") {
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
    # group F - S. anginosus, incl. S. anginosus anginosus and S. anginosus whileyi
    out[out %like_case% "^B_STRPT_ANGN(_|$)"] <- "B_STRPT_GRPF"
    # group G - S. dysgalactiae and S. canis (though dysgalactiae is also group C and will be matched there)
    out[out %like_case% "^B_STRPT_(DYSG|CANS)(_|$)"] <- "B_STRPT_GRPG"
    # group H - S. sanguinis
    out[out %like_case% "^B_STRPT_SNGN(_|$)"] <- "B_STRPT_GRPH"
    # group K - S. salivarius, incl. S. salivarius salivariuss and S. salivarius thermophilus
    out[out %like_case% "^B_STRPT_SLVR(_|$)"] <- "B_STRPT_GRPK"
    # group L - only S. dysgalactiae which is also group C & G, so ignore it here
  }

  # All unknowns ----
  out[is.na(out) & !is.na(x)] <- "UNKNOWN"
  AMR_env$mo_failures <- unique(x[out == "UNKNOWN" & !x %in% c("UNKNOWN", "con") & !x %like_case% "^[(]unknown [a-z]+[)]$" & !is.na(x)])
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
  paste0(
    "(",
    "[^A-Za-z- \\(\\)\\[\\]{}]+",
    "|",
    "([({]|\\[).+([})]|\\])",
    "|",
    "(^| )(e?spp|e?ssp|e?ss|e?sp|e?subsp|sube?species|biovar|biotype|serovar|serogr.?up|e?species)[.]*( |$|(complex|group)$))"
  )
}

# UNDOCUMENTED METHODS ----------------------------------------------------

# will be exported using s3_register() in R/zzz.R
pillar_shaft.mo <- function(x, ...) {
  add_MO_lookup_to_AMR_env()
  out <- format(x)
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

# will be exported using s3_register() in R/zzz.R
type_sum.mo <- function(x, ...) {
  "mo"
}

# will be exported using s3_register() in R/zzz.R
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
          big.mark = ",",
          decimal.mark = "."
        ),
        " (", percentage(sum(grams == "Gram-negative", na.rm = TRUE) / length(grams),
          digits = digits
        ),
        ")"
      ),
      `Gram-positive` = paste0(
        format(sum(grams == "Gram-positive", na.rm = TRUE),
          big.mark = ",",
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

# will be exported using s3_register() in R/zzz.R
get_skimmers.mo <- function(column) {
  skimr::sfl(
    skim_type = "mo",
    unique_total = ~ length(unique(stats::na.omit(.))),
    gram_negative = ~ sum(mo_is_gram_negative(.), na.rm = TRUE),
    gram_positive = ~ sum(mo_is_gram_positive(.), na.rm = TRUE),
    top_genus = ~ names(sort(-table(mo_genus(stats::na.omit(.), language = NULL))))[1L],
    top_species = ~ names(sort(-table(mo_name(stats::na.omit(.), language = NULL))))[1L]
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
print.mo_uncertainties <- function(x, ...) {
  if (NROW(x) == 0) {
    cat(word_wrap("No uncertainties to show. Only uncertainties of the last call of `as.mo()` or any `mo_*()` function are stored.\n\n", add_fn = font_blue))
    return(invisible(NULL))
  }

  cat(word_wrap("Matching scores are based on the resemblance between the input and the full taxonomic name, and the pathogenicity in humans. See `?mo_matching_score`.\n\n", add_fn = font_blue))
  
  add_MO_lookup_to_AMR_env()
  
  col_red <- function(x) font_rose_bg(font_black(x, collapse = NULL), collapse = NULL)
  col_orange <- function(x) font_orange_bg(font_black(x, collapse = NULL), collapse = NULL)
  col_yellow <- function(x) font_yellow_bg(font_black(x, collapse = NULL), collapse = NULL)
  col_green <- function(x) font_green_bg(font_black(x, collapse = NULL), collapse = NULL)
  
  if (has_colour()) {
    cat(word_wrap("Colour keys: ",
      col_red(" 0.000-0.499 "),
      col_orange(" 0.500-0.599 "),
      col_yellow(" 0.600-0.699 "),
      col_green(" 0.700-1.000"),
      add_fn = font_blue
    ), font_green_bg(" "), "\n", sep = "")
  }
  
  score_set_colour <- function(text, scores) {
    # set colours to scores
    text[scores >= 0.7] <- col_green(text[scores >= 0.7])
    text[scores >= 0.6 & scores < 0.7] <- col_yellow(text[scores >= 0.6 & scores < 0.7])
    text[scores >= 0.5 & scores < 0.6] <- col_orange(text[scores >= 0.5 & scores < 0.6])
    text[scores < 0.5] <- col_red(text[scores < 0.5])
    text
  }

  txt <- ""
  for (i in seq_len(nrow(x))) {
    if (x[i, ]$candidates != "") {
      candidates <- unlist(strsplit(x[i, ]$candidates, ", ", fixed = TRUE))
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
          ),
          ifelse(n_candidates == 25,
            font_grey(" [showing first 25]"),
            ""
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
      # Add "Based on {input}" text if it differs from the original input
      ifelse(x[i, ]$original_input != x[i, ]$input, paste0(strrep(" ", nchar(x[i, ]$original_input) + 6), "Based on input \"", x[i, ]$input, "\""), ""),
      # Add note if result was coerced to accepted taxonomic name
      ifelse(x[i, ]$keep_synonyms == FALSE & x[i, ]$mo %in% AMR_env$MO_lookup$mo[which(AMR_env$MO_lookup$status == "synonym")],
        paste0(
          strrep(" ", nchar(x[i, ]$original_input) + 6),
          font_red(paste0("This old taxonomic name was converted to ", font_italic(AMR_env$MO_lookup$fullname[match(synonym_mo_to_accepted_mo(x[i, ]$mo), AMR_env$MO_lookup$mo)], collapse = NULL), " (", synonym_mo_to_accepted_mo(x[i, ]$mo), ")."), collapse = NULL)
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
  out[x %like_case% "^g[abcdfghkl]s$"] <- gsub("g([abcdfghkl])s",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "^g[abcdfghkl]s$"],
    perl = TRUE
  )
  # Streptococci in different languages, like "estreptococos grupo B"
  out[x %like_case% "strepto[ck]o[ck].* [abcdfghkl]$"] <- gsub(".*e?strepto[ck]o[ck].* ([abcdfghkl])$",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "strepto[ck]o[ck].* [abcdfghkl]$"],
    perl = TRUE
  )
  out[x %like_case% "strep[a-z]* group [abcdfghkl]$"] <- gsub(".* ([abcdfghkl])$",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "strep[a-z]* group [abcdfghkl]$"],
    perl = TRUE
  )
  out[x %like_case% "group [abcdfghkl] strepto[ck]o[ck]"] <- gsub(".*group ([abcdfghkl]) strepto[ck]o[ck].*",
    "B_STRPT_GRP\\U\\1",
    x[x %like_case% "group [abcdfghkl] strepto[ck]o[ck]"],
    perl = TRUE
  )
  out[x %like_case% "ha?emoly.*strep"] <- "B_STRPT_HAEM"
  out[x %like_case% "(strepto.* mil+er+i|^mgs[^a-z]*$)"] <- "B_STRPT_MILL"
  out[x %like_case% "mil+er+i gr"] <- "B_STRPT_MILL"
  out[x %like_case% "((strepto|^s).* viridans|^vgs[^a-z]*$)"] <- "B_STRPT_VIRI"
  out[x %like_case% "(viridans.* (strepto|^s).*|^vgs[^a-z]*$)"] <- "B_STRPT_VIRI"

  # Salmonella in different languages, like "Salmonella grupo B"
  out[x %like_case% "salmonella.* [bcd]$"] <- gsub(".*salmonella.* ([bcd])$",
    "B_SLMNL_GRP\\U\\1",
    x[x %like_case% "salmonella.* [bcd]$"],
    perl = TRUE
  )
  out[x %like_case% "group [bcd] salmonella"] <- gsub(".*group ([bcd]) salmonella*",
    "B_SLMNL_GRP\\U\\1",
    x[x %like_case% "group [bcd] salmonella"],
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
  out[is.na(out) & x %like_case% "anaerob[a-z]+ (micro)?.*organism"] <- "B_ANAER"

  # yeasts and fungi
  out[x %like_case% "^yeast?"] <- "F_YEAST"
  out[x %like_case% "^fung(us|i)"] <- "F_FUNGUS"

  # trivial names known to the field
  out[x %like_case% "meningo[ck]o[ck]"] <- "B_NESSR_MNNG"
  out[x %like_case% "gono[ck]o[ck]"] <- "B_NESSR_GNRR"
  out[x %like_case% "pneumo[ck]o[ck]"] <- "B_STRPT_PNMN"

  # unexisting names (xxx and con are WHONET codes)
  out[x %in% c("con", "other", "none", "unknown") | x %like_case% "virus"] <- "UNKNOWN"

  # WHONET has a lot of E. coli and Vibrio cholerae names
  out[x %like_case% "escherichia coli"] <- "B_ESCHR_COLI"
  out[x %like_case% "vibrio cholerae"] <- "B_VIBRI_CHLR"

  out
}

italicise <- function(x) {
  out <- font_italic(x, collapse = NULL)
  out[x %like_case% "Salmonella [A-Z]"] <- paste(
    font_italic("Salmonella"),
    gsub("Salmonella ", "", x[x %like_case% "Salmonella [A-Z]"])
  )
  out[x %like_case% "Streptococcus [A-Z]"] <- paste(
    font_italic("Streptococcus"),
    gsub("Streptococcus ", "", x[x %like_case% "Streptococcus [A-Z]"])
  )
  if (has_colour()) {
    out <- gsub("(Group|group|Complex|complex)(\033\\[23m)?", "\033[23m\\1", out, perl = TRUE)
  }
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
      error = function(e) stop(e$message, call. = FALSE)
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
  reference_df[, "mo"] <- as.mo(reference_df[, "mo", drop = TRUE])
  reference_df
}

get_mo_uncertainties <- function() {
  remember <- list(uncertainties = AMR_env$mo_uncertainties)
  # empty them, otherwise e.g. mo_shortname("Chlamydophila psittaci") will give 3 notes
  AMR_env$mo_uncertainties <- NULL
  remember
}

load_mo_uncertainties <- function(metadata) {
  AMR_env$mo_uncertainties <- metadata$uncertainties
}

synonym_mo_to_accepted_mo <- function(x, fill_in_accepted = FALSE) {
  x_gbif <- AMR_env$MO_lookup$gbif_renamed_to[match(x, AMR_env$MO_lookup$mo)]
  x_lpsn <- AMR_env$MO_lookup$lpsn_renamed_to[match(x, AMR_env$MO_lookup$mo)]
  x_gbif[!x_gbif %in% AMR_env$MO_lookup$gbif] <- NA
  x_lpsn[!x_lpsn %in% AMR_env$MO_lookup$lpsn] <- NA

  out <- ifelse(is.na(x_lpsn),
    AMR_env$MO_lookup$mo[match(x_gbif, AMR_env$MO_lookup$gbif)],
    AMR_env$MO_lookup$mo[match(x_lpsn, AMR_env$MO_lookup$lpsn)]
  )
  if (isTRUE(fill_in_accepted)) {
    x_accepted <- which(AMR_env$MO_lookup$status[match(x, AMR_env$MO_lookup$mo)] == "accepted")
    out[x_accepted] <- x[x_accepted]
  }
  out
}

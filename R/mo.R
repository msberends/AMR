# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

#' Transform Input to a Microorganism Code
#'
#' Use this function to determine a valid microorganism code ([`mo`]). Determination is done using intelligent rules and the complete taxonomic kingdoms `r vector_and(unique(microorganisms$kingdom[which(!grepl("(unknown|Fungi)", microorganisms$kingdom))]), quotes = FALSE)`, and most microbial species from the kingdom Fungi (see *Source*). The input can be almost anything: a full name (like `"Staphylococcus aureus"`), an abbreviated name (such as `"S. aureus"`), an abbreviation known in the field (such as `"MRSA"`), or just a genus. See *Examples*.
#' @param x a [character] vector or a [data.frame] with one or two columns
#' @param Becker a [logical] to indicate whether staphylococci should be categorised into coagulase-negative staphylococci ("CoNS") and coagulase-positive staphylococci ("CoPS") instead of their own species, according to Karsten Becker *et al.* (1,2,3).
#'
#' This excludes *Staphylococcus aureus* at default, use `Becker = "all"` to also categorise *S. aureus* as "CoPS".
#' @param Lancefield a [logical] to indicate whether a beta-haemolytic *Streptococcus* should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield (4). These streptococci will be categorised in their first group, e.g. *Streptococcus dysgalactiae* will be group C, although officially it was also categorised into groups G and L.
#'
#' This excludes enterococci at default (who are in group D), use `Lancefield = "all"` to also categorise all enterococci as group D.
#' @param minimum_matching_score a numeric value to set as the lower limit for the [MO matching score][mo_matching_score()]. When left blank, this will be determined automatically based on the character length of `x`, its [taxonomic kingdom][microorganisms] and [human pathogenicity][mo_matching_score()].
#' @param allow_uncertain a number between `0` (or `"none"`) and `3` (or `"all"`), or `TRUE` (= `2`) or `FALSE` (= `0`) to indicate whether the input should be checked for less probable results, see *Details*
#' @param reference_df a [data.frame] to be used for extra reference when translating `x` to a valid [`mo`]. See [set_mo_source()] and [get_mo_source()] to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ignore_pattern a regular expression (case-insensitive) of which all matches in `x` must return `NA`. This can be convenient to exclude known non-relevant input and can also be set with the option `AMR_ignore_pattern`, e.g. `options(AMR_ignore_pattern = "(not reported|contaminated flora)")`.
#' @param language language to translate text like "no growth", which defaults to the system language (see [get_AMR_locale()])
#' @param info a [logical] to indicate if a progress bar should be printed if more than 25 items are to be coerced, defaults to `TRUE` only in interactive mode
#' @param ... other arguments passed on to functions
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
#' ## General Info
#'
#' A microorganism (MO) code from this package (class: [`mo`]) is human readable and typically looks like these examples:
#' ```
#'   Code               Full name
#'   ---------------    --------------------------------------
#'   B_KLBSL            Klebsiella
#'   B_KLBSL_PNMN       Klebsiella pneumoniae
#'   B_KLBSL_PNMN_RHNS  Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |    |
#'   |   |    |    |
#'   |   |    |    \---> subspecies, a 4-5 letter acronym
#'   |   |    \----> species, a 4-5 letter acronym
#'   |   \----> genus, a 5-7 letter acronym
#'   \----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                             F (Fungi), PL (Plantae), P (Protozoa)
#' ```
#'
#' Values that cannot be coerced will be considered 'unknown' and will get the MO code `UNKNOWN`.
#'
#' Use the [`mo_*`][mo_property()] functions to get properties based on the returned code, see *Examples*.
#'
#' The algorithm uses data from the List of Prokaryotic names with Standing in Nomenclature (LPSN) and the Global Biodiversity Information Facility (GBIF) (see [microorganisms]).
#'
#' The [as.mo()] function uses several coercion rules for fast and logical results. It assesses the input matching criteria in the following order:
#'
#' 1. Human pathogenic prevalence: the function  starts with more prevalent microorganisms, followed by less prevalent ones;
#' 2. Taxonomic kingdom: the function starts with determining Bacteria, then Fungi, then Protozoa, then others;
#' 3. Breakdown of input values to identify possible matches.
#'
#' This will lead to the effect that e.g. `"E. coli"` (a microorganism highly prevalent in humans) will return the microbial ID of *Escherichia coli* and not *Entamoeba coli* (a microorganism less prevalent in humans), although the latter would alphabetically come first.
#'
#' ## Coping with Uncertain Results
#'
#' In addition, the [as.mo()] function can differentiate four levels of uncertainty to guess valid results:
#' - Uncertainty level 0: no additional rules are applied;
#' - Uncertainty level 1: allow previously accepted (but now invalid) taxonomic names and minor spelling errors;
#' - Uncertainty level 2: allow all of level 1, strip values between brackets, inverse the words of the input, strip off text elements from the end keeping at least two elements;
#' - Uncertainty level 3: allow all of level 1 and 2, strip off text elements from the end, allow any part of a taxonomic name.
#'
#' The level of uncertainty can be set using the argument `allow_uncertain`. The default is `allow_uncertain = TRUE`, which is equal to uncertainty level 2. Using `allow_uncertain = FALSE` is equal to uncertainty level 0 and will skip all rules. You can also use e.g. `as.mo(..., allow_uncertain = 1)` to only allow up to level 1 uncertainty.
#'
#' With the default setting (`allow_uncertain = TRUE`, level 2), below examples will lead to valid results:
#' - `"Streptococcus group B (known as S. agalactiae)"`. The text between brackets will be removed and a warning will be thrown that the result *Streptococcus group B* (``r as.mo("Streptococcus group B")``) needs review.
#' - `"S. aureus - please mind: MRSA"`. The last word will be stripped, after which the function will try to find a match. If it does not, the second last word will be stripped, etc. Again, a warning will be thrown that the result *Staphylococcus aureus* (``r as.mo("Staphylococcus aureus")``) needs review.
#' - `"Fluoroquinolone-resistant Neisseria gonorrhoeae"`. The first word will be stripped, after which the function will try to find a match. A warning will be thrown that the result *Neisseria gonorrhoeae* (``r as.mo("Neisseria gonorrhoeae")``) needs review.
#'
#' There are three helper functions that can be run after using the [as.mo()] function:
#' - Use [mo_uncertainties()] to get a [data.frame] that prints in a pretty format with all taxonomic names that were guessed. The output contains the matching score for all matches (see *Matching Score for Microorganisms* below).
#' - Use [mo_failures()] to get a [character] [vector] with all values that could not be coerced to a valid value.
#' - Use [mo_renamed()] to get a [data.frame] with all values that could be coerced based on old, previously accepted taxonomic names.
#'
#' ## Microbial Prevalence of Pathogens in Humans
#'
#' The coercion rules consider the prevalence of microorganisms in humans grouped into three groups, which is available as the `prevalence` columns in the [microorganisms] data set. The grouping into human pathogenic prevalence is explained in the section *Matching Score for Microorganisms* below.
#' @inheritSection mo_matching_score Matching Score for Microorganisms
#  (source as a section here, so it can be inherited by other man pages:)
#' @section Source:
#' 1. Becker K. *et al.* (2014). **Coagulase-Negative Staphylococci.** *Clin Microbiol Rev.* 27(4): 870-926; \doi{10.1128/CMR.00109-13}
#' 2. Becker K. *et al.* (2019). **Implications of identifying the recently defined members of the *S. aureus* complex, *S. argenteus* and *S. schweitzeri*: A position paper of members of the ESCMID Study Group for staphylococci and Staphylococcal Diseases (ESGS).** *Clin Microbiol Infect*; \doi{10.1016/j.cmi.2019.02.028}
#' 3. Becker K. *et al.* (2020). **Emergence of coagulase-negative staphylococci** *Expert Rev Anti Infect Ther.* 18(4):349-366; \doi{10.1080/14787210.2020.1730813}
#' 4. Lancefield R.C. (1933). **A serological differentiation of human and other groups of hemolytic streptococci**. *J Exp Med.* 57(4): 571-95; \doi{10.1084/jem.57.4.571}
#' 5. Berends M.S. *et al.* (2022). **Trends in Occurrence and Phenotypic Resistance of Coagulase-Negative Staphylococci (CoNS) Found in Human Blood in the Northern Netherlands between 2013 and 2019** *Microorganisms* 10(9), 1801; \doi{10.3390/microorganisms10091801}
#' 6. `r TAXONOMY_VERSION$LPSN$citation` Accessed from <`r TAXONOMY_VERSION$LPSN$url`> on `r documentation_date(TAXONOMY_VERSION$LPSN$accessed_date)`.
#' 7. `r TAXONOMY_VERSION$GBIF$citation` Accessed from <`r TAXONOMY_VERSION$GBIF$url`> on `r documentation_date(TAXONOMY_VERSION$GBIF$accessed_date)`.
#' 8. `r TAXONOMY_VERSION$SNOMED$citation` URL: <`r TAXONOMY_VERSION$SNOMED$url`>
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
#'   "Staphylococcus aureus",
#'   "Staphylococcus aureus (MRSA,",
#'   "Zthafilokkoockus oureuz", # handles incorrect spelling
#'   "MRSA", # Methicillin Resistant S. aureus
#'   "VISA", # Vancomycin Intermediate S. aureus
#'   "VRSA", # Vancomycin Resistant S. aureus
#'   115329001
#' )) # SNOMED CT code
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
#' mo_genus("Esch coli")
#' mo_gramstain("E. coli")
#' mo_is_intrinsic_resistant("E. coli", "vanco")
#' }
as.mo <- function(x,
                  Becker = FALSE,
                  Lancefield = FALSE,
                  minimum_matching_score = NULL,
                  allow_uncertain = TRUE,
                  keep_synonyms = FALSE,
                  reference_df = get_mo_source(),
                  ignore_pattern = getOption("AMR_ignore_pattern"),
                  language = get_AMR_locale(),
                  info = interactive(),
                  ...) {
  meet_criteria(x, allow_class = c("mo", "data.frame", "list", "character", "numeric", "integer", "factor"), allow_NA = TRUE)
  meet_criteria(Becker, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(Lancefield, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(keep_synonyms, allow_class = c("logical", "character"), has_length = 1)
  meet_criteria(minimum_matching_score, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE)
  meet_criteria(reference_df, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(ignore_pattern, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  language <- validate_language(language)
  meet_criteria(info, allow_class = "logical", has_length = 1)

  if (tryCatch(all(x[!is.na(x)] %in% AMR::microorganisms$mo) &
    isFALSE(Becker) &
    isTRUE(keep_synonyms) &&
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

  # WHONET: xxx = no growth
  x[tolower(x) %in% c("", "xxx", "na", "nan")] <- NA_character_

  if (tryCatch(all(x == "" | gsub(".*(unknown ).*", "unknown name", tolower(x), perl = TRUE) %in% MO_lookup$fullname_lower, na.rm = TRUE) &&
    isFALSE(Becker) &&
    isTRUE(keep_synonyms) &&
    isFALSE(Lancefield), error = function(e) FALSE)) {
    # to improve speed, special case for taxonomically correct full names (case-insensitive)
    return(set_clean_class(MO_lookup[match(
      gsub(".*(unknown ).*", "unknown name",
        tolower(x),
        perl = TRUE
      ),
      MO_lookup$fullname_lower
    ), "mo", drop = TRUE],
    new_class = c("mo", "character")
    ))
  }

  out <- rep(NA_character_, length(x))

  # below we use base R's match(), known for powering '%in%', and incredibly fast!

  # From reference_df ----
  reference_df <- repair_reference_df(reference_df)
  if (!is.null(reference_df)) {
    out[x %in% reference_df[[1]]] <- reference_df[[2]][match(x[x %in% reference_df[[1]]], reference_df[[1]])]
  }
  # From MO code ----
  out[is.na(out) & x %in% AMR::microorganisms$mo] <- x[is.na(out) & x %in% AMR::microorganisms$mo]
  # From full name ----
  out[is.na(out) & x %in% AMR::microorganisms$fullname] <- AMR::microorganisms$mo[match(x[is.na(out) & x %in% AMR::microorganisms$fullname], AMR::microorganisms$fullname)]
  # From known codes ----
  out[is.na(out) & x %in% AMR::microorganisms.codes$code] <- AMR::microorganisms.codes$mo[match(x[is.na(out) & x %in% AMR::microorganisms.codes$code], AMR::microorganisms.codes$code)]
  # From SNOMED ----
  if (any(is.na(out) & x %in% unlist(microorganisms$snomed), na.rm = TRUE)) {
    # found this extremely fast gem here: https://stackoverflow.com/a/11002456/4575331
    out[is.na(out) & x %in% unlist(microorganisms$snomed)] <- microorganisms$mo[rep(seq_along(microorganisms$snomed), vapply(FUN.VALUE = double(1), microorganisms$snomed, length))[match(x[is.na(out) & x %in% unlist(microorganisms$snomed)], unlist(microorganisms$snomed))]]
  }
  # From previous hits in this session ----
  old <- out
  out[is.na(out) & x %in% pkg_env$mo_previously_coerced$x] <- pkg_env$mo_previously_coerced$mo[match(x[is.na(out) & x %in% pkg_env$mo_previously_coerced$x], pkg_env$mo_previously_coerced$x)]
  new <- out
  if (isTRUE(info) && message_not_thrown_before("as.mo", old[seq_len(min(100, length(old)))], new[seq_len(min(100, length(new)))], entire_session = TRUE) && any(is.na(old) & !is.na(new), na.rm = TRUE)) {
    message_(
      "Returning previously coerced value", ifelse(sum(is.na(old) & !is.na(new)) > 1, "s", ""),
      " for ", vector_and(x[is.na(old) & !is.na(new)]), ". Run `mo_reset_session()` to reset this."
    )
  }

  # For all other input ----
  if (any(is.na(out) & !is.na(x))) {
    # reset uncertainties
    pkg_env$mo_uncertainties <- pkg_env$mo_uncertainties[0, ]

    # Laboratory systems: remove (translated) entries like "no growth", "not E. coli", etc.
    x[trimws2(x) %like% translate_into_language("no .*growth", language = language)] <- NA_character_
    x[trimws2(x) %like% paste0("^(", translate_into_language("no|not", language = language), ") ")] <- NA_character_

    # run over all unique leftovers
    x_unique <- unique(x[is.na(out) & !is.na(x)])

    # set up progress bar
    progress <- progress_ticker(n = length(x_unique), n_min = 10, print = info)
    on.exit(close(progress))

    # run it
    x_coerced <- lapply(x_unique, function(x_search) {
      progress$tick()

      x_out <- trimws(tolower(x_search))
      x_parts <- strsplit(gsub("-", " ", x_out, fixed = TRUE), " ", fixed = TRUE)[[1]]

      # do a pre-match on first character (and if it contains a space, first chars of first two terms)
      if (length(x_parts) == 2) {
        filtr <- which(MO_lookup$full_first == substr(x_parts[1], 1, 1) & MO_lookup$species_first == substr(x_parts[2], 1, 1))
      } else if (length(x_parts) > 2) {
        first_chars <- paste0("(^| )", "[", paste(substr(x_parts, 1, 1), collapse = ""), "]")
        filtr <- which(MO_lookup$full_first %like_case% first_chars)
      } else if (nchar(x_out) == 4) {
        # no space and 4 characters - probably a code such as STAU or ESCO!
        if (isTRUE(info)) {
          message_("Input \"", x_search, "\" is assumed to be a microorganism code - trying to match on ", vector_and(c(substr(x_out, 1, 2), substr(x_out, 3, 4)), sort = FALSE))
        }
        filtr <- which(MO_lookup$fullname_lower %like_case% paste0("(^| )", substr(x_out, 1, 2), ".* ", substr(x_out, 3, 4)))
      } else if (nchar(x_out) <= 6) {
        # no space and 5-6 characters - probably a code such as STAAUR or ESCCOL!
        first_part <- paste0(substr(x_out, 1, 2), "[a-z]*", substr(x_out, 3, 3))
        second_part <- substr(x_out, 4, nchar(x_out))
        if (isTRUE(info)) {
          message_("Input \"", x_search, "\" is assumed to be a microorganism code - trying to match on ", vector_and(c(gsub("[a-z]*", "(...)", first_part, fixed = TRUE), second_part), sort = FALSE))
        }
        filtr <- which(MO_lookup$fullname_lower %like_case% paste0("(^| )", first_part, ".* ", second_part))
      } else {
        filtr <- which(MO_lookup$full_first == substr(x_out, 1, 1))
      }
      if (length(filtr) == 0) {
        mo_to_search <- MO_lookup$fullname
      } else {
        mo_to_search <- MO_lookup$fullname[filtr]
      }
      pkg_env$mo_to_search <- mo_to_search
      # determine the matching score on the original search value
      m <- mo_matching_score(x = x_search, n = mo_to_search)
      if (is.null(minimum_matching_score)) {
        minimum_matching_score_current <- min(0.7, min(10, nchar(x_search)) * 0.08)
        # correct back for prevalence
        minimum_matching_score_current <- minimum_matching_score_current / MO_lookup$prevalence[match(mo_to_search, MO_lookup$fullname)]
        # correct back for kingdom
        minimum_matching_score_current <- minimum_matching_score_current / MO_lookup$kingdom_index[match(mo_to_search, MO_lookup$fullname)]
      } else {
        minimum_matching_score_current <- minimum_matching_score
      }
      m[m < minimum_matching_score_current] <- NA_real_
      top_hits <- mo_to_search[order(m, decreasing = TRUE, na.last = NA)] # na.last = NA will remove the NAs
      if (length(top_hits) == 0) {
        warning_("No hits found for \"", x_search, "\" with minimum_matching_score = ", ifelse(is.null(minimum_matching_score), "NULL", minimum_matching_score), ". Try setting this value higher.")
        result_mo <- NA_character_
      } else {
        result_mo <- MO_lookup$mo[match(top_hits[1], MO_lookup$fullname)]
        pkg_env$mo_uncertainties <- rbind(pkg_env$mo_uncertainties,
          data.frame(
            minimum_matching_score = ifelse(is.null(minimum_matching_score), "NULL", minimum_matching_score),
            input = x_search,
            fullname = top_hits[1],
            mo = result_mo,
            candidates = ifelse(length(top_hits) > 1, paste(top_hits[2:min(26, length(top_hits))], collapse = ", "), ""),
            stringsAsFactors = FALSE
          ),
          stringsAsFactors = FALSE
        )
        # save to package env to save time for next time
        pkg_env$mo_previously_coerced <- unique(rbind(pkg_env$mo_previously_coerced,
          data.frame(
            x = x_search,
            mo = result_mo,
            stringsAsFactors = FALSE
          ),
          stringsAsFactors = FALSE
        ))
      }
      # the actual result:
      result_mo
    })

    # remove progress bar from console
    close(progress)
    # expand from unique again
    out[is.na(out)] <- unlist(x_coerced)[match(x[is.na(out)], x_unique)]

    # Throw note about uncertainties ----
    if (isTRUE(info) && NROW(pkg_env$mo_uncertainties) > 0) {
      if (message_not_thrown_before("as.mo", "uncertainties", pkg_env$mo_uncertainties$input)) {
        plural <- c("", "this")
        if (length(pkg_env$mo_uncertainties$input) > 1) {
          plural <- c("s", "these uncertainties")
        }
        if (length(pkg_env$mo_uncertainties$input) <= 3) {
          examples <- vector_and(paste0(
            '"', pkg_env$mo_uncertainties$input,
            '" (assumed ', font_italic(pkg_env$mo_uncertainties$fullname, collapse = NULL), ")"
          ),
          quotes = FALSE
          )
        } else {
          examples <- paste0(nr2char(length(pkg_env$mo_uncertainties$input)), " microorganism", plural[1])
        }
        msg <- paste0(
          "Microorganism translation was uncertain for ", examples,
          ". Run `mo_uncertainties()` to review ", plural[2], "."
        )
        message_(msg)
      }
    }
  } # end of loop over all yet unknowns

  # Keep or replace synonyms ----
  if (isFALSE(keep_synonyms)) {
    out_old <- out

    gbif_matches <- AMR::microorganisms$gbif_renamed_to[match(out, AMR::microorganisms$mo)]
    gbif_matches[!gbif_matches %in% AMR::microorganisms$gbif] <- NA
    out[which(!is.na(gbif_matches))] <- AMR::microorganisms$mo[match(gbif_matches[which(!is.na(gbif_matches))], AMR::microorganisms$gbif)]

    lpsn_matches <- AMR::microorganisms$lpsn_renamed_to[match(out, AMR::microorganisms$mo)]
    lpsn_matches[!lpsn_matches %in% AMR::microorganisms$lpsn] <- NA
    out[which(!is.na(lpsn_matches))] <- AMR::microorganisms$mo[match(lpsn_matches[which(!is.na(lpsn_matches))], AMR::microorganisms$lpsn)]

    if (isTRUE(info) && (any(!is.na(gbif_matches)) || any(!is.na(lpsn_matches))) && message_not_thrown_before("as.mo", gbif_matches[which(!is.na(gbif_matches))][1:5], lpsn_matches[which(!is.na(lpsn_matches))][1:5]) && length(c(lpsn_matches, gbif_matches)) > 0) {
      total_old <- out_old[which(!is.na(gbif_matches) | !is.na(lpsn_matches))]
      total_new <- out[which(!is.na(gbif_matches) | !is.na(lpsn_matches))]

      total_new <- total_new[!duplicated(total_old)]
      total_old <- total_old[!duplicated(total_old)]

      total_new <- total_new[order(total_old)]
      total_old <- total_old[order(total_old)]

      refs_old <- microorganisms$ref[match(total_old, microorganisms$mo)]
      refs_old[!is.na(refs_old)] <- paste0(" (", refs_old[!is.na(refs_old)], ")")
      refs_old[is.na(refs_old)] <- ""
      refs_new <- microorganisms$ref[match(total_new, microorganisms$mo)]
      refs_new[!is.na(refs_new)] <- paste0(" (", refs_new[!is.na(refs_new)], ")")
      refs_new[is.na(refs_new)] <- ""

      message_(
        "The following microorganism", ifelse(length(total_old) > 1, "s were", " was"), " taxonomically renamed (use `keep_synonyms = TRUE` to leave uncorrected):\n",
        paste0("  ", microorganisms$fullname[match(total_old, microorganisms$mo)],
          refs_old,
          " -> ", microorganisms$fullname[match(total_new, microorganisms$mo)],
          refs_new,
          collapse = "\n"
        )
      )
    }
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
    if (any(out %in% AMR::microorganisms$mo[match(post_Becker, AMR::microorganisms$fullname)])) {
      if (message_not_thrown_before("as.mo", "becker")) {
        warning_("in `as.mo()`: Becker ", font_italic("et al."), " (2014, 2019, 2020) does not contain these species named after their publication: ",
          vector_and(font_italic(gsub("Staphylococcus", "S.", post_Becker, fixed = TRUE), collapse = NULL), quotes = FALSE),
          ". Categorisation to CoNS/CoPS was taken from the original scientific publication(s).",
          immediate = TRUE
        )
      }
    }

    # 'MO_CONS' and 'MO_COPS' are <mo> vectors created in R/_pre_commit_hook.R
    out[out %in% MO_CONS] <- "B_STPHY_CONS"
    out[out %in% MO_COPS] <- "B_STPHY_COPS"
    if (Becker == "all") {
      out[out == "B_STPHY_AURS"] <- "B_STPHY_COPS"
    }
  }

  # Apply Lancefield ----
  if (isTRUE(Lancefield) || Lancefield == "all") {
    # group A - S. pyogenes
    out[out == "B_STRPT_PYGN"] <- "B_STRPT_GRPA"
    # group B - S. agalactiae
    out[out == "B_STRPT_AGLC"] <- "B_STRPT_GRPB"
    # group C - all subspecies within S. dysgalactiae and S. equi (such as S. equi zooepidemicus)
    out[out %like_case% "^B_STRPT_(DYSG|EQUI)(_|$)"] <- "B_STRPT_GRPC"
    if (Lancefield == "all") {
      # group D - all enterococci
      out[out %like_case% "^B_ENTRC(_|$)"] <- "B_STRPT_GRPD"
    }
    # group F - S. anginosus, incl. S. anginosus anginosus and S. anginosus whileyi
    out[out %like_case% "^B_STRPT_ANGN(_|$)"] <- "B_STRPT_GRPF"
    # group G - only S. dysgalactiae which is also group C, so ignore it here
    # group H - S. sanguinis
    out[out == "B_STRPT_SNGN"] <- "B_STRPT_GRPH"
    # group K - S. salivarius, incl. S. salivarius salivariuss and S. salivarius thermophilus
    out[out %like_case% "^B_STRPT_SLVR(_|$)"] <- "B_STRPT_GRPK"
    # group L - only S. dysgalactiae which is also group C, so ignore it here
  }

  # Return class ----
  set_clean_class(out,
    new_class = c("mo", "character")
  )
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  inherits(x, "mo")
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.mo <- function(x, ...) {
  out <- format(x)
  # grey out the kingdom (part until first "_")
  out[!is.na(x)] <- gsub("^([A-Z]+_)(.*)", paste0(font_subtle("\\1"), "\\2"), out[!is.na(x)], perl = TRUE)
  # and grey out every _
  out[!is.na(x)] <- gsub("_", font_subtle("_"), out[!is.na(x)])

  # markup NA and UNKNOWN
  out[is.na(x)] <- font_na("  NA")
  out[x == "UNKNOWN"] <- font_na("  UNKNOWN")

  df <- tryCatch(get_current_data(arg_name = "x", call = 0),
    error = function(e) NULL
  )
  if (!is.null(df)) {
    mo_cols <- vapply(FUN.VALUE = logical(1), df, is.mo)
  } else {
    mo_cols <- NULL
  }

  if (!all(x[!is.na(x)] %in% MO_lookup$mo) |
    (!is.null(df) && !all(unlist(df[, which(mo_cols), drop = FALSE]) %in% MO_lookup$mo))) {
    # markup old mo codes
    out[!x %in% MO_lookup$mo] <- font_italic(font_na(x[!x %in% MO_lookup$mo],
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
      "Please update your MO codes with `as.mo()`."
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
  cat("Class <mo>\n")
  x_names <- names(x)
  if (is.null(x_names) & print.shortnames == TRUE) {
    x_names <- tryCatch(mo_shortname(x, ...), error = function(e) NULL)
  }
  x <- as.character(x)
  names(x) <- x_names
  if (!all(x[!is.na(x)] %in% MO_lookup$mo)) {
    warning_(
      "Some MO codes are from a previous AMR package version. ",
      "Please update the MO codes with `as.mo()`."
    )
  }
  print.default(x, quote = FALSE)
}

#' @method summary mo
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- as.mo(object) # force again, could be mo from older pkg version
  top <- as.data.frame(table(x), responseName = "n", stringsAsFactors = FALSE)
  top_3 <- top[order(-top$n), 1, drop = TRUE][1:3]
  value <- c(
    "Class" = "mo",
    "<NA>" = length(x[is.na(x)]),
    "Unique" = pm_n_distinct(x[!is.na(x)]),
    "#1" = top_3[1],
    "#2" = top_3[2],
    "#3" = top_3[3]
  )
  class(value) <- c("summaryDefault", "table")
  value
}

#' @method as.data.frame mo
#' @export
#' @noRd
as.data.frame.mo <- function(x, ...) {
  if (!all(x[!is.na(x)] %in% MO_lookup$mo)) {
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
  return_after_integrity_check(y, "microorganism code", as.character(AMR::microorganisms$mo))
}
#' @method [[<- mo
#' @export
#' @noRd
"[[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  # must only contain valid MOs
  return_after_integrity_check(y, "microorganism code", as.character(AMR::microorganisms$mo))
}
#' @method c mo
#' @export
#' @noRd
c.mo <- function(...) {
  x <- list(...)[[1L]]
  y <- NextMethod()
  attributes(y) <- attributes(x)
  return_after_integrity_check(y, "microorganism code", as.character(AMR::microorganisms$mo))
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

#' @rdname as.mo
#' @export
mo_failures <- function() {
  pkg_env$mo_failures
}

#' @rdname as.mo
#' @export
mo_uncertainties <- function() {
  if (is.null(pkg_env$mo_uncertainties)) {
    return(NULL)
  }
  set_clean_class(as.data.frame(pkg_env$mo_uncertainties,
    stringsAsFactors = FALSE
  ),
  new_class = c("mo_uncertainties", "data.frame")
  )
}

#' @method print mo_uncertainties
#' @export
#' @noRd
print.mo_uncertainties <- function(x, ...) {
  if (NROW(x) == 0) {
    return(NULL)
  }
  cat(word_wrap("Matching scores are based on pathogenicity in humans and the resemblance between the input and the full taxonomic name. See `?mo_matching_score`.\n\n", add_fn = font_blue))
  if (has_colour()) {
    cat(word_wrap("Colour keys: ",
      font_red_bg(" 0.000-0.499 "),
      font_orange_bg(" 0.500-0.599 "),
      font_yellow_bg(" 0.600-0.699 "),
      font_green_bg(" 0.700-1.000"),
      add_fn = font_blue
    ), font_green_bg(" "), "\n", sep = "")
  }

  score_set_colour <- function(text, scores) {
    # set colours to scores
    text[scores >= 0.7] <- font_green_bg(text[scores >= 0.7], collapse = NULL)
    text[scores >= 0.6 & scores < 0.7] <- font_yellow_bg(text[scores >= 0.6 & scores < 0.7], collapse = NULL)
    text[scores >= 0.5 & scores < 0.6] <- font_orange_bg(text[scores >= 0.5 & scores < 0.6], collapse = NULL)
    text[scores < 0.5] <- font_red_bg(text[scores < 0.5], collapse = NULL)
    text
  }

  txt <- ""
  for (i in seq_len(nrow(x))) {
    if (x[i, ]$candidates != "") {
      candidates <- unlist(strsplit(x[i, ]$candidates, ", ", fixed = TRUE))
      scores <- mo_matching_score(x = x[i, ]$input, n = candidates)
      n_candidates <- length(candidates)

      candidates_formatted <- font_italic(candidates, collapse = NULL)
      scores_formatted <- trimws(formatC(round(scores, 3), format = "f", digits = 3))
      scores_formatted <- score_set_colour(scores_formatted, scores)

      # sort on descending scores
      candidates_formatted <- candidates_formatted[order(1 - scores)]
      scores_formatted <- scores_formatted[order(1 - scores)]

      candidates <- word_wrap(paste0(
        "Also matched: ",
        vector_and(paste0(
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
      extra_indent = nchar("Also matched: ")
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
        strwrap(
          paste0(
            '"', x[i, ]$input, '"',
            " -> ",
            paste0(
              font_bold(font_italic(x[i, ]$fullname)),
              ifelse(!is.na(x[i, ]$renamed_to), paste(", renamed to", font_italic(x[i, ]$renamed_to)), ""),
              " (", x[i, ]$mo,
              ", ", score_set_colour(score_formatted, score),
              ") "
            )
          ),
          width = 0.98 * getOption("width"),
          exdent = nchar(x[i, ]$input) + 6
        ),
        collapse = "\n"
      ),
      candidates,
      sep = "\n"
    )
    txt <- paste0(gsub("\n\n", "\n", txt), "\n\n")
  }
  cat(txt)
}

#' @rdname as.mo
#' @export
mo_reset_session <- function() {
  if (NROW(pkg_env$mo_previously_coerced) > 0) {
    message_("Reset ", NROW(pkg_env$mo_previously_coerced), " previously matched input values.")
    pkg_env$mo_previously_coerced <- pkg_env$mo_previously_coerced[0, , drop = FALSE]
  } else {
    message_("No previously matched input values to reset.")
  }
}

#' @rdname as.mo
#' @export
mo_renamed <- function() {
  items <- pkg_env$mo_renamed
  if (is.null(items)) {
    items <- data.frame(stringsAsFactors = FALSE)
  } else {
    items <- pm_distinct(items, old_name, .keep_all = TRUE)
  }
  set_clean_class(as.data.frame(items,
    stringsAsFactors = FALSE
  ),
  new_class = c("mo_renamed", "data.frame")
  )
}

#' @method print mo_renamed
#' @export
#' @noRd
print.mo_renamed <- function(x, ...) {
  if (NROW(x) == 0) {
    return(invisible())
  }
  for (i in seq_len(nrow(x))) {
    message_(
      font_italic(x$old_name[i]),
      ifelse(x$old_ref[i] %in% c("", NA),
        "",
        paste0(" (", gsub("et al.", font_italic("et al."), x$old_ref[i]), ")")
      ),
      " was renamed ",
      ifelse(!x$new_ref[i] %in% c("", NA) && as.integer(gsub("[^0-9]", "", x$new_ref[i])) < as.integer(gsub("[^0-9]", "", x$old_ref[i])),
        font_bold("back to "),
        ""
      ),
      font_italic(x$new_name[i]),
      ifelse(x$new_ref[i] %in% c("", NA),
        "",
        paste0(" (", gsub("et al.", font_italic("et al."), x$new_ref[i]), ")")
      ),
      " [", x$mo[i], "]"
    )
  }
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
    stop_ifnot(allow_uncertain %in% c(0:3),
      '`allow_uncertain` must be a number between 0 (or "none") and 3 (or "all"), or TRUE (= 2) or FALSE (= 0)',
      call = FALSE
    )
  }
  allow_uncertain
}

get_mo_failures_uncertainties_renamed <- function() {
  remember <- list(
    failures = pkg_env$mo_failures,
    uncertainties = pkg_env$mo_uncertainties,
    renamed = pkg_env$mo_renamed
  )
  # empty them, otherwise mo_shortname("Chlamydophila psittaci") will give 3 notes
  pkg_env$mo_failures <- NULL
  pkg_env$mo_uncertainties <- NULL
  pkg_env$mo_renamed <- NULL
  remember
}

load_mo_failures_uncertainties_renamed <- function(metadata) {
  pkg_env$mo_failures <- metadata$failures
  pkg_env$mo_uncertainties <- metadata$uncertainties
  pkg_env$mo_renamed <- metadata$renamed
}

trimws2 <- function(x) {
  trimws(gsub("[\\s]+", " ", x, perl = TRUE))
}

parse_and_convert <- function(x) {
  tryCatch(
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
      parsed <- gsub(" +", " ", parsed, perl = TRUE)
      parsed <- trimws(parsed)
      parsed
    },
    error = function(e) stop(e$message, call. = FALSE)
  ) # this will also be thrown when running `as.mo(no_existing_object)`
  parsed
}

replace_old_mo_codes <- function(x, property) {
  # this function transform old MO codes to current codes, such as:
  # B_ESCH_COL (AMR v0.5.0) -> B_ESCHR_COLI
  ind <- x %like_case% "^[A-Z]_[A-Z_]+$" & !x %in% MO_lookup$mo
  if (any(ind)) {
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
        results <- MO_lookup$mo[MO_lookup$kingdom %like_case% kingdom &
          MO_lookup$fullname_lower %like_case% name]
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

strip_words <- function(text, n, side = "right") {
  out <- lapply(strsplit(text, " "), function(x) {
    if (side %like% "^r" & length(x) > n) {
      x[seq_len(length(x) - n)]
    } else if (side %like% "^l" & length(x) > n) {
      x[2:length(x)]
    }
  })
  vapply(FUN.VALUE = character(1), out, paste, collapse = " ")
}

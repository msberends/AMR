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
# Center Groningen in the Netherlands, in collaboration with many      #
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


# SESSION CACHES ------------------------------------------------------
#
# as.ab() and as.mo() share one cache design, kept in the package environment `AMR_env`:
#
# 1. Lookup indices, derived from `AMR_env$AB_lookup` and `AMR_env$MO_lookup` (e.g. all ATC codes or SNOMED
#    codes as one flat vector). These are built lazily on first use, see get_ab_index() and get_mo_index().
# 2. Coercion caches, `AMR_env$ab_previously_coerced` and `AMR_env$mo_previously_coerced`, that remember the
#    outcome of earlier coercions. Both have the same structure: a `key` (the cleaned input together with the
#    settings that affect the outcome), the original `input`, and the resulting code as `value`.
# 3. Translated names of antimicrobials and microorganisms, built once per language, see get_ab_translated_index() and
#    get_mo_translated_index().
# 4. For microorganisms, the mapping from (outdated) codes to codes of currently accepted names, which is
#    memoised per unique code, see synonym_mo_to_accepted_mo().
#
# All of these are invalidated with reset_ab_cache() or reset_mo_cache(), which must be called whenever the
# lookup tables change, e.g. in add_custom_antimicrobials() and add_custom_microorganisms().
# The prepared translation tables of translate_into_language(), see get_translation_table(), are kept for the whole
# session, as the internal translations never change.

new_coercion_cache <- function() {
  data.frame(
    key = character(0),
    input = character(0),
    value = character(0),
    stringsAsFactors = FALSE
  )
}

coercion_cache_key <- function(x, ...) {
  # all arguments in `...` are settings that can change the outcome of a coercion, they become part of the key
  if (length(x) == 0) {
    return(character(0))
  }
  settings <- vapply(FUN.VALUE = character(1), list(...), function(s) paste0(s, collapse = ""))
  paste(x, paste(settings, collapse = "|"), sep = "|")
}

coercion_cache_has <- function(type, keys) {
  keys %in% AMR_env[[paste0(type, "_previously_coerced")]]$key
}

coercion_cache_get <- function(type, keys) {
  cache <- AMR_env[[paste0(type, "_previously_coerced")]]
  cache$value[match(keys, cache$key)]
}

coercion_cache_add <- function(type, keys, inputs, values) {
  cache_name <- paste0(type, "_previously_coerced")
  keep <- !is.na(keys) & !duplicated(keys)
  if (!any(keep)) {
    return(invisible(NULL))
  }
  keys <- keys[keep]
  cache <- AMR_env[[cache_name]]
  AMR_env[[cache_name]] <- rbind_AMR(
    cache[which(!cache$key %in% keys), , drop = FALSE],
    data.frame(
      key = keys,
      input = as.character(inputs[keep]),
      value = as.character(values[keep]),
      stringsAsFactors = FALSE
    )
  )
  invisible(NULL)
}

reset_ab_cache <- function() {
  AMR_env$AB_index <- NULL
  AMR_env$AB_translated <- NULL
  AMR_env$custom_ab_id_keys <- NULL
  AMR_env$ab_previously_coerced <- new_coercion_cache()
}

mo_codes_valid <- function(x) {
  # whether codes are in `AMR_env$MO_lookup`; a hashed set avoids scanning the whole lookup table for every check, as `[<-.mo` runs one for every piece that rbind() binds
  x <- as.character(x)
  if (length(x) > 1000) {
    # for long vectors, base matching is faster than looking up every unique value
    add_MO_lookup_to_AMR_env()
    return(x %in% as.character(AMR_env$MO_lookup$mo))
  }
  if (is.null(AMR_env$MO_codes_hash)) {
    add_MO_lookup_to_AMR_env()
    codes <- unique(as.character(AMR_env$MO_lookup$mo))
    AMR_env$MO_codes_hash <- list2env(stats::setNames(as.list(rep(TRUE, length(codes))), codes), hash = TRUE, size = length(codes))
  }
  ux <- unique(x)
  valid <- vapply(
    FUN.VALUE = logical(1),
    ux,
    function(code) !is.na(code) && nzchar(code) && exists(code, envir = AMR_env$MO_codes_hash, inherits = FALSE),
    USE.NAMES = FALSE
  )
  valid[match(x, ux)]
}

reset_mo_cache <- function() {
  AMR_env$MO_codes_hash <- NULL
  AMR_env$MO_index <- NULL
  AMR_env$MO_translated <- NULL
  AMR_env$mo_accepted <- NULL
  AMR_env$mo_previously_coerced <- new_coercion_cache()
  AMR_env$mo_previously_uncertain <- NULL
}

get_ab_index <- function() {
  if (is.null(AMR_env$AB_index)) {
    ab <- as.character(AMR_env$AB_lookup$ab)
    atc <- unlist(AMR_env$AB_lookup$atc, use.names = FALSE)
    atc_ab <- rep(ab, lengths(AMR_env$AB_lookup$atc))
    synonyms <- tolower(unlist(AMR_env$AB_lookup$synonyms, use.names = FALSE))
    synonyms_ab <- rep(ab, lengths(AMR_env$AB_lookup$synonyms))
    # generalised LOINC codes, synonyms, and abbreviations as flat vectors, with the row number in `AMR_env$AB_lookup`
    flat <- function(col) {
      vals <- AMR_env$AB_lookup[[col]]
      list(value = unlist(vals, use.names = FALSE), row = rep(seq_along(vals), lengths(vals)))
    }
    gen_syn <- flat("generalised_synonyms")
    # synonyms of at least 5 characters, for the Levenshtein distance in as.ab()
    gen_syn_long <- which(!is.na(gen_syn$value) & nchar(gen_syn$value) >= 5)
    AMR_env$AB_index <- list(
      atc = atc[!is.na(atc)],
      atc_ab = atc_ab[!is.na(atc)],
      synonyms = synonyms[!is.na(synonyms)],
      synonyms_ab = synonyms_ab[!is.na(synonyms)],
      generalised_loinc = flat("generalised_loinc"),
      generalised_synonyms = gen_syn,
      generalised_synonyms_long = list(value = gen_syn$value[gen_syn_long], row = gen_syn$row[gen_syn_long]),
      generalised_abbreviations = flat("generalised_abbreviations")
    )
  }
  AMR_env$AB_index
}

get_ab_translated_index <- function(language) {
  # translated and generalised names of all antimicrobials, used by as.ab() for input that could not be matched otherwise
  # this is built once per language
  if (is.null(AMR_env$AB_translated[[language]])) {
    AMR_env$AB_translated[[language]] <- generalise_antibiotic_name(translate_AMR(AMR_env$AB_lookup$name, language = language))
  }
  AMR_env$AB_translated[[language]]
}

get_mo_index <- function() {
  add_MO_lookup_to_AMR_env()
  if (is.null(AMR_env$MO_index)) {
    snomed <- unlist(AMR_env$MO_lookup$snomed, use.names = FALSE)
    snomed_mo <- rep(as.character(AMR_env$MO_lookup$mo), lengths(AMR_env$MO_lookup$snomed))
    AMR_env$MO_index <- list(
      snomed = snomed[!is.na(snomed)],
      snomed_mo = snomed_mo[!is.na(snomed)]
    )
  }
  AMR_env$MO_index
}

get_mo_translated_index <- function(language) {
  # translated names of microorganisms, such as "Staphylocoque à coagulase négative (CoNS)", so that as.mo() recognises the output of mo_name() in that language
  # this is built once per language and only contains the few names that change by translation
  add_MO_lookup_to_AMR_env()
  if (is.null(AMR_env$MO_translated[[language]])) {
    trns <- TRANSLATIONS[which(TRANSLATIONS$affect_mo_name == TRUE), , drop = FALSE]
    fullname <- AMR_env$MO_lookup$fullname
    is_translatable <- grepl(paste0("(", paste(trns$pattern, collapse = "|"), ")"), fullname, perl = TRUE)
    fullname <- fullname[is_translatable]
    mo <- as.character(AMR_env$MO_lookup$mo[is_translatable])
    translated <- translate_into_language(fullname, language = language, only_affect_mo_names = TRUE)
    is_changed <- translated != fullname
    translated <- translated[is_changed]
    mo <- mo[is_changed]
    # also without the abbreviation, such as "(CoNS)"
    translated_short <- gsub(" +\\([^()]+\\)$", "", translated, perl = TRUE)
    name <- tolower(c(translated, translated_short))
    mo <- c(mo, mo)
    AMR_env$MO_translated[[language]] <- list(
      name = name[!duplicated(name)],
      mo = mo[!duplicated(name)]
    )
  }
  AMR_env$MO_translated[[language]]
}

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
# 3. For microorganisms, the mapping from (outdated) codes to codes of currently accepted names, which is
#    memoised per unique code, see synonym_mo_to_accepted_mo().
#
# All of these are invalidated with reset_ab_cache() or reset_mo_cache(), which must be called whenever the
# lookup tables change, e.g. in add_custom_antimicrobials() and add_custom_microorganisms().

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
  AMR_env$custom_ab_id_keys <- NULL
  AMR_env$ab_previously_coerced <- new_coercion_cache()
}

reset_mo_cache <- function() {
  AMR_env$MO_index <- NULL
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
    AMR_env$AB_index <- list(
      atc = atc[!is.na(atc)],
      atc_ab = atc_ab[!is.na(atc)],
      synonyms = synonyms[!is.na(synonyms)],
      synonyms_ab = synonyms_ab[!is.na(synonyms)]
    )
  }
  AMR_env$AB_index
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

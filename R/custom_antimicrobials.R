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

#' Add Custom Antimicrobials
#'
#' With [add_custom_antimicrobials()] you can add your own custom antimicrobial drug names and codes.
#' @param x A [data.frame] resembling the [antimicrobials] data set, at least containing columns "ab" and "name". To add synonyms to existing antimicrobials, columns "ab" and "synonyms" suffice (see *Details*).
#' @param ab For [add_custom_antimicrobial_synonyms()]: a [character] vector of existing antimicrobial codes (see [as.ab()]), or a [data.frame] with columns "ab" and "synonyms" (one synonym per row).
#' @param synonyms A [character] vector of synonyms. If `ab` has length 1, all `synonyms` are added to that antimicrobial; otherwise `ab` and `synonyms` must have the same length.
#' @details **Important:** Due to how \R works, the [add_custom_antimicrobials()] function has to be run in every \R session - added antimicrobials are not stored between sessions and are thus lost when \R is exited.
#'
#' There are two ways to circumvent this and automate the process of adding antimicrobials:
#'
#' **Method 1:** Using the package option [`AMR_custom_ab`][AMR-options], which is the preferred method. To use this method:
#'
#'    1. Create a data set in the structure of the [antimicrobials] data set (containing at the very least columns "ab" and "name", or only columns "ab" and "synonyms" for synonyms of existing antimicrobials) and save it with [saveRDS()] to a location of choice, e.g. `"~/my_custom_ab.rds"`, or any remote location.
#'
#'    2. Set the file location to the package option [`AMR_custom_ab`][AMR-options]: `options(AMR_custom_ab = "~/my_custom_ab.rds")`. This can even be a remote file location, such as an https URL. Since options are not saved between \R sessions, it is best to save this option to the `.Rprofile` file so that it will be loaded on start-up of \R. To do this, open the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`, add this text and save the file:
#'
#'       ```r
#'       # Add custom antimicrobial codes:
#'       options(AMR_custom_ab = "~/my_custom_ab.rds")
#'       ```
#'
#'       Upon package load, this file will be loaded and run through the [add_custom_antimicrobials()] function.
#'
#'       This file can also hold synonyms for antimicrobials that already exist, such as local trade names: rows with an existing code in column "ab" and only column "synonyms" filled (one synonym per row, codes may repeat) are passed on to [add_custom_antimicrobial_synonyms()]. Rows with an existing code that have any other column filled still give an error. Synonyms must be valid UTF-8 text; convert other encodings first, e.g. with `iconv(x, from = "CP949", to = "UTF-8")`.
#'
#' **Method 2:** Loading the antimicrobial additions directly from your `.Rprofile` file. Note that the definitions will be stored in a user-specific \R file, which is a suboptimal workflow. To use this method:
#'
#'    1. Edit the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`.
#'
#'    2. Add a text like below and save the file:
#'
#'       ```r
#'        # Add custom antibiotic drug codes:
#'        AMR::add_custom_antimicrobials(
#'          data.frame(ab = "TESTAB",
#'                     name = "Test Antibiotic",
#'                     group = "Test Group")
#'        )
#'       ```
#'
#' Use [add_custom_antimicrobial_synonyms()] to add extra names, such as local trade names, to antimicrobials that already exist (including antimicrobials added with [add_custom_antimicrobials()]). These synonyms are recognised by [as.ab()] and all `ab_*()` functions. Synonyms may contain non-ASCII characters (e.g. trade names in other scripts): they are matched against the input as given, before [as.ab()] transliterates the input to ASCII. Matching ignores case, all white space (including no-break and ideographic spaces) and zero-width characters. If the input does not match as given, a trailing part in parentheses or brackets (such as `(piperacillin, tazobactam)`) and a trailing strength written as a number with a unit symbol (such as `4.5g`, `4/0.5 g` or `500 mg/mL`) are removed one at a time, and matching is tried again after each removal, so that the longest registered synonym wins. Full-width characters are read as their ASCII equivalents. No words are removed, in any language: to recognise a name with a dosage form (such as `Inj`), add that variant as a synonym as well. Only user-added synonyms are matched this way, there is no other fuzzy matching for them, and a strength on its own is never matched to a user-added synonym. A synonym cannot be added if it already identifies another antimicrobial as its code, name, synonym, abbreviation, ATC code, CID or LOINC code, also when written with other punctuation (such as `J01-CR05`).
#'
#' Use [clear_custom_antimicrobials()] to clear the previously added antimicrobials and synonyms.
#' @seealso [add_custom_microorganisms()] to add custom microorganisms.
#' @rdname add_custom_antimicrobials
#' @export
#' @examples
#' \donttest{
#' # returns a wildly guessed result:
#' as.ab("testab")
#'
#' # now add a custom entry - it will be considered by as.ab() and
#' # all ab_*() functions
#' add_custom_antimicrobials(
#'   data.frame(
#'     ab = "TESTAB",
#'     name = "Test Antibiotic",
#'     # you can add any property present in the
#'     # 'antimicrobials' data set, such as 'group':
#'     group = "Test Group"
#'   )
#' )
#'
#' # "testab" is now a new antibiotic:
#' as.ab("testab")
#' ab_name("testab")
#' ab_group("testab")
#'
#' ab_info("testab")
#'
#'
#' # Add Co-fluampicil, which is one of the many J01CR50 codes, see
#' # https://atcddd.fhi.no/ddd/list_of_ddds_combined_products/
#' add_custom_antimicrobials(
#'   data.frame(
#'     ab = "COFLU",
#'     name = "Co-fluampicil",
#'     atc = "J01CR50",
#'     group = "Beta-lactams/penicillins"
#'   )
#' )
#' ab_atc("Co-fluampicil")
#' ab_name("J01CR50")
#'
#' # even antimicrobial selectors work
#' # see ?amr_selector
#' x <- data.frame(
#'   random_column = "some value",
#'   coflu = as.sir("S"),
#'   ampicillin = as.sir("R")
#' )
#' x
#' x[, betalactams()]
#'
#' # add local trade names as synonyms to existing antimicrobials,
#' # here a Korean trade name of piperacillin/tazobactam
#' add_custom_antimicrobial_synonyms("TZP", "\uD0C0\uC870\uC2E0\uC8FC")
#' ab_name("\uD0C0\uC870\uC2E0\uC8FC")
#' # a trailing strength, as in many hospital exports, is ignored
#' as.ab("\uD0C0\uC870\uC2E0\uC8FC 4.5g")
#' }
add_custom_antimicrobials <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  # remove any extra class/type, such as grouped tbl, or data.table:
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  stop_ifnot(
    "ab" %in% colnames(x),
    "`x` must contain columns \"ab\" and \"name\", or columns \"ab\" and \"synonyms\" to add synonyms to existing antimicrobials."
  )
  # rows that only fill in "ab" and "synonyms" add synonyms to existing antimicrobials, see add_custom_antimicrobial_synonyms()
  synonym_rows <- custom_ab_synonym_rows(x)
  records <- x[!synonym_rows, , drop = FALSE]
  stop_if(
    NROW(records) > 0 && !"name" %in% colnames(records),
    "`x` must contain columns \"ab\" and \"name\", or columns \"ab\" and \"synonyms\" to add synonyms to existing antimicrobials."
  )
  stop_if(
    any(records$ab %in% AMR_env$AB_lookup$ab),
    "Antimicrobial drug code(s) ", vector_and(records$ab[records$ab %in% AMR_env$AB_lookup$ab]), " already exist in the internal `antimicrobials` data set."
  )
  # names of new antimicrobials must not already be in use as user-added synonyms of other antimicrobials,
  # compared in the same ways as custom_ab_synonyms_prepare() does, so that the order of adding does not matter
  if (NROW(records) > 0 && NROW(AMR_env$custom_ab_synonyms) > 0) {
    syn_tab <- AMR_env$custom_ab_synonyms
    rec_ids <- custom_ab_identifier_keys(records)
    rec_gen <- custom_ab_generalised_punct(generalise_antibiotic_name(as.character(records$name)))
    hits <- data.frame(
      rec = c(rec_ids$ab, rec_ids$ab, as.character(records$ab)),
      s = c(
        match(rec_ids$key, syn_tab$key, incomparables = NA),
        match(rec_ids$pkey, custom_ab_punct_key(syn_tab$key), incomparables = NA),
        match(rec_gen, custom_ab_generalised_punct(custom_ab_generalised(syn_tab$synonym)), incomparables = NA)
      ),
      stringsAsFactors = FALSE
    )
    hits <- hits[!is.na(hits$s), , drop = FALSE]
    clash <- syn_tab$ab[hits$s] != hits$rec
    stop_if(
      any(clash),
      "The name or synonym(s) ", vector_and(unique(syn_tab$synonym[hits$s[clash]])), " of the new antimicrobial(s) are already in use as user-added synonyms of ", vector_and(unique(syn_tab$ab[hits$s[clash]])), "."
    )
  }
  # validate the synonyms before anything is added, so that an error leaves the session unchanged
  syn <- NULL
  if (any(synonym_rows)) {
    syn_col <- x$synonyms[synonym_rows]
    syn_ab <- as.character(x$ab[synonym_rows])
    if (is.list(syn_col)) {
      syn_ab <- rep(syn_ab, vapply(FUN.VALUE = integer(1), syn_col, length))
      syn_col <- unlist(syn_col, use.names = FALSE)
    }
    syn <- custom_ab_synonyms_prepare(ab = syn_ab, synonyms = syn_col, new_records = records)
  }
  if (NROW(records) > 0) {
    add_custom_antimicrobials_records(records)
  }
  if (!is.null(syn)) {
    custom_ab_synonyms_commit(syn)
  }
  invisible(NULL)
}

# rows of `x` that have a "synonyms" column and no other column filled besides "ab"
custom_ab_synonym_rows <- function(x) {
  if (!"synonyms" %in% colnames(x)) {
    return(rep(FALSE, NROW(x)))
  }
  other_filled <- rep(FALSE, NROW(x))
  for (col in setdiff(colnames(x), c("ab", "synonyms"))) {
    other_filled <- other_filled | custom_ab_is_filled(x[, col, drop = TRUE])
  }
  !other_filled
}

custom_ab_is_filled <- function(col) {
  if (is.list(col)) {
    vapply(FUN.VALUE = logical(1), col, function(v) any(!is.na(v) & as.character(v) != ""))
  } else {
    !is.na(col) & as.character(col) != ""
  }
}

add_custom_antimicrobials_records <- function(x) {
  # keep only columns available in the antimicrobials data set
  x <- x[, colnames(AMR_env$AB_lookup)[colnames(AMR_env$AB_lookup) %in% colnames(x)], drop = FALSE]
  x$generalised_name <- generalise_antibiotic_name(x$name)
  x$generalised_all <- as.list(x$generalised_name)
  for (col in colnames(x)) {
    if (is.list(AMR_env$AB_lookup[, col, drop = TRUE]) & !is.list(x[, col, drop = TRUE])) {
      x[, col] <- as.list(x[, col, drop = TRUE])
    }
  }

  AMR_env$custom_ab_codes <- c(AMR_env$custom_ab_codes, x$ab)
  class(AMR_env$AB_lookup$ab) <- "character"

  new_df <- AMR_env$AB_lookup[0, , drop = FALSE][seq_len(NROW(x)), , drop = FALSE]
  rownames(new_df) <- NULL
  list_cols <- vapply(FUN.VALUE = logical(1), new_df, is.list)
  for (l in which(list_cols)) {
    # prevent binding NULLs in lists, replace with NA
    new_df[, l] <- as.list(NA_character_)
  }
  for (col in colnames(x)) {
    # assign new values
    new_df[, col] <- x[, col, drop = TRUE]
  }
  AMR_env$AB_lookup <- unique(rbind_AMR(AMR_env$AB_lookup, new_df))

  class(AMR_env$AB_lookup$ab) <- c("ab", "character")
  # the lookup table changed, so earlier coercions and lookup indices are no longer valid
  reset_ab_cache()
  message_("Added ", nr2char(nrow(x)), " record", ifelse(nrow(x) > 1, "s", ""), " to the internal {.code antimicrobials} data set.")
}

#' @rdname add_custom_antimicrobials
#' @export
add_custom_antimicrobial_synonyms <- function(ab, synonyms = NULL) {
  if (is.data.frame(ab)) {
    stop_ifnot(
      all(c("ab", "synonyms") %in% colnames(ab)),
      "`ab` must contain columns \"ab\" and \"synonyms\" when it is a data.frame."
    )
    synonyms <- ab$synonyms
    ab <- as.character(ab$ab)
    if (is.list(synonyms)) {
      ab <- rep(ab, vapply(FUN.VALUE = integer(1), synonyms, length))
      synonyms <- unlist(synonyms, use.names = FALSE)
    }
  } else {
    meet_criteria(ab, allow_class = c("character", "ab"))
    meet_criteria(synonyms, allow_class = "character")
    ab <- as.character(ab)
    if (length(ab) == 1) {
      ab <- rep(ab, length(synonyms))
    }
    stop_ifnot(
      length(ab) == length(synonyms),
      "`ab` must be of length 1 or of the same length as `synonyms`."
    )
  }
  syn <- custom_ab_synonyms_prepare(ab = ab, synonyms = synonyms)
  if (is.null(syn)) {
    message_("No synonyms to add.")
    return(invisible(NULL))
  }
  custom_ab_synonyms_commit(syn)
  invisible(NULL)
}

# checks synonyms and returns a data.frame(ab, synonym, key), or NULL if there is nothing to add;
# `new_records` are antimicrobials that are about to be added in the same call of add_custom_antimicrobials()
custom_ab_synonyms_prepare <- function(ab, synonyms, new_records = NULL) {
  ab <- as.character(ab)
  synonyms <- enc2utf8(as.character(synonyms))
  invalid <- !is.na(synonyms) & !custom_ab_valid_utf8(synonyms)
  stop_if(
    any(invalid),
    "Synonym(s) must be valid UTF-8 text, which does not apply to the synonym(s) at position ", vector_and(which(invalid), quotes = FALSE), ". Convert them first, e.g. with `iconv(x, from = \"CP949\", to = \"UTF-8\")`."
  )
  synonyms <- trimws2(synonyms)
  keep <- !is.na(ab) & !is.na(synonyms) & !is.na(custom_ab_synonym_key(synonyms))
  ab <- ab[keep]
  synonyms <- synonyms[keep]
  if (length(ab) == 0) {
    return(NULL)
  }
  new_codes <- if (is.null(new_records)) character(0) else as.character(new_records$ab)
  known <- c(as.character(AMR_env$AB_lookup$ab), new_codes)
  unknown <- unique(ab[!ab %in% known])
  stop_if(
    length(unknown) > 0,
    "Antimicrobial drug code(s) ", vector_and(unknown), " do not exist. Use add_custom_antimicrobials() to add new antimicrobials first."
  )
  keys <- custom_ab_synonym_key(synonyms)
  # a synonym must not already identify another antimicrobial - its code, name, synonyms, abbreviations,
  # ATC codes, CID or LOINC codes - not with other punctuation (e.g. "J01-CR05") and not in the
  # generalised form that as.ab() uses for its exact matches either
  ids <- custom_ab_identifier_keys_cached()
  gen_ids <- custom_ab_generalised_ids()
  if (!is.null(new_records) && NROW(new_records) > 0) {
    ids <- rbind(ids, custom_ab_identifier_keys(new_records))
    new_gen <- data.frame(ab = as.character(new_records$ab), gen = generalise_antibiotic_name(as.character(new_records$name)), stringsAsFactors = FALSE)
    new_gen$genp <- custom_ab_generalised_punct(new_gen$gen)
    gen_ids <- rbind(gen_ids, new_gen[!is.na(new_gen$genp), , drop = FALSE])
  }
  owners_of <- function(table_keys, table_ab, k) {
    groups <- split(table_ab, table_keys)
    lapply(match(k, names(groups), incomparables = NA), function(j) if (is.na(j)) character(0) else groups[[j]])
  }
  by_key <- owners_of(ids$key, ids$ab, keys)
  by_pkey <- owners_of(ids$pkey, ids$ab, custom_ab_punct_key(keys))
  by_gen <- owners_of(gen_ids$genp, gen_ids$ab, custom_ab_generalised_punct(custom_ab_generalised(synonyms)))
  owners <- lapply(seq_along(keys), function(i) setdiff(unique(c(by_key[[i]], by_pkey[[i]], by_gen[[i]])), ab[i]))
  taken <- lengths(owners) > 0
  stop_if(
    any(taken),
    "Synonym(s) already identify another antimicrobial (as its code, name, synonym, abbreviation, ATC code, CID or LOINC code): ",
    vector_and(unique(paste0(synonyms[taken], " (", vapply(FUN.VALUE = character(1), owners[taken], paste, collapse = ", "), ")")), quotes = FALSE), "."
  )
  existing <- AMR_env$custom_ab_synonyms
  clash <- keys %in% existing$key & ab != existing$ab[match(keys, existing$key)]
  clash <- clash | (duplicated(keys) & !duplicated(paste(keys, ab)))
  stop_if(
    any(clash),
    "Synonym(s) ", vector_and(unique(synonyms[clash])), " would refer to more than one antimicrobial."
  )
  data.frame(ab = ab, synonym = synonyms, key = keys, stringsAsFactors = FALSE)
}

custom_ab_synonyms_commit <- function(syn) {
  class(AMR_env$AB_lookup$ab) <- "character"
  for (code in unique(syn$ab)) {
    i <- which(AMR_env$AB_lookup$ab == code)[1L]
    new_syn <- syn$synonym[syn$ab == code]
    old_syn <- AMR_env$AB_lookup$synonyms[[i]]
    old_syn <- old_syn[!is.na(old_syn) & old_syn != ""]
    # only added to `synonyms` (for ab_synonyms()); matching is done by the exact match at the start of as.ab(),
    # so these synonyms are deliberately kept out of the fuzzy matching on `generalised_*`
    AMR_env$AB_lookup$synonyms[[i]] <- unique(c(old_syn, new_syn))
  }
  class(AMR_env$AB_lookup$ab) <- c("ab", "character")
  AMR_env$custom_ab_synonyms <- unique(rbind(AMR_env$custom_ab_synonyms, syn))
  # the lookup table changed, so earlier coercions and lookup indices are no longer valid
  reset_ab_cache()
  n_syn <- length(unique(syn$synonym))
  n_ab <- length(unique(syn$ab))
  message_(
    "Added ", nr2char(n_syn), " synonym", ifelse(n_syn > 1, "s", ""),
    " to ", nr2char(n_ab), " antimicrobial", ifelse(n_ab > 1, "s", ""), "."
  )
}

# validUTF8() requires R >= 3.3.0, so use iconv(), which returns NA for invalid input
custom_ab_valid_utf8 <- function(x) {
  !is.na(iconv(x, from = "UTF-8", to = "UTF-8"))
}

# white space for the user-added synonym functions: all Unicode white space, like trimws2()
custom_ab_ws <- "[\\h\\v\\p{Z}\u200B\u200C\u200D\u2060\uFEFF]"

# key for user-added synonyms: case and all white space are ignored, full-width characters are read as ASCII,
# other characters are kept as given;
# input that is not valid UTF-8 or that consists of white space only gets an NA key and thus never matches
custom_ab_synonym_key <- function(x) {
  x <- enc2utf8(as.character(x))
  out <- rep(NA_character_, length(x))
  valid <- !is.na(x) & custom_ab_valid_utf8(x)
  out[valid] <- toupper(gsub(paste0(custom_ab_ws, "+"), "", custom_ab_fullwidth_to_ascii(x[valid]), perl = TRUE))
  out[!is.na(out) & out == ""] <- NA_character_
  out
}

# generalised form as used by the exact stage of as.ab(), or NA when the text cannot be transliterated to ASCII
custom_ab_generalised <- function(x) {
  out <- rep(NA_character_, length(x))
  ascii <- suppressWarnings(iconv(toupper(enc2utf8(as.character(x))), from = "UTF-8", to = "ASCII//TRANSLIT"))
  ok <- !is.na(ascii) & !grepl("?", ascii, fixed = TRUE)
  if (any(ok)) {
    # some iconv implementations transliterate diacritics to separate marks (e.g. "O" with acute to "'O")
    gen <- generalise_antibiotic_name(gsub("['`^~\"]", "", ascii[ok]))
    gen[!grepl("[A-Z0-9]", gen)] <- NA_character_
    out[ok] <- gen
  }
  out
}

# key without punctuation (letters and digits only), so that a synonym cannot take over another
# antimicrobial's identifier by adding punctuation, e.g. "J01-CR05" or "Zo-syn"
custom_ab_punct_key <- function(key) {
  out <- gsub("[^\\p{L}\\p{N}]+", "", key, perl = TRUE)
  out[!is.na(out) & out == ""] <- NA_character_
  out
}

# the same for generalised forms, which consist of A-Z, 0-9 and punctuation only
custom_ab_generalised_punct <- function(gen) {
  out <- gsub("[^A-Z0-9]+", "", gen)
  out[!is.na(out) & out == ""] <- NA_character_
  out
}

# parts after a product name that are removed one at a time by custom_ab_strip_step(): trailing punctuation,
# a part in parentheses or brackets (also nested) and a strength written as a number with a unit symbol;
# this is deliberately language-independent, so no words (such as dosage forms) are removed
custom_ab_strip_patterns <- local({
  num <- "[0-9]+([.,][0-9]+)*"
  unit <- "(mg|g|mcg|ug|\u00B5g|\u03BCg|ng|kg|iu|miu|u|ml|l|meq|mmol|%|\u338E|\u338D|\u338F|\u3396|\u2113)"
  sep <- "[ ,_*-]*"
  per <- paste0("( ?/ ?(", num, ")? ?", unit, ")?")
  combination <- paste0("(", num, " ?", unit, "? ?[/:] ?)?")
  c(
    punctuation = "[.,;:]+$",
    parenthesised = " ?(\\((?:[^()]++|(?1))*\\))$",
    bracketed = " ?(\\[(?:[^][]++|(?1))*\\])$",
    strength = paste0("(?<![0-9.,])", sep, combination, num, " ?", unit, per, "\\.?$")
  )
})

# full-width ASCII characters to ASCII
custom_ab_fullwidth_to_ascii <- function(x) {
  wide <- !is.na(x) & grepl("[\uFF01-\uFF5E]", x, perl = TRUE)
  if (any(wide)) {
    x[wide] <- vapply(FUN.VALUE = character(1), x[wide], function(s) {
      cp <- utf8ToInt(s)
      w <- cp >= 65281L & cp <= 65374L
      cp[w] <- cp[w] - 65248L
      intToUtf8(cp)
    }, USE.NAMES = FALSE)
  }
  x
}

# white space to single spaces and full-width ASCII characters to ASCII
custom_ab_strip_normalise <- function(x) {
  x <- custom_ab_fullwidth_to_ascii(enc2utf8(as.character(x)))
  trimws2(gsub(paste0(custom_ab_ws, "+"), " ", x, perl = TRUE))
}

# removes one trailing part (see custom_ab_strip_patterns) from each value, e.g. "Brand (comment) 4.5g" ->
# "Brand (comment)" -> "Brand"; values that cannot be shortened (or would become empty) are returned as is;
# only used by as.ab() to look up user-added synonyms when the input as given does not match one
custom_ab_strip_step <- function(x) {
  out <- x
  todo <- !is.na(x)
  for (p in custom_ab_strip_patterns) {
    if (!any(todo)) break
    y <- trimws2(gsub(p, "", out[todo], perl = TRUE, ignore.case = TRUE))
    changed <- y != out[todo] & y != ""
    if (any(changed)) {
      i <- which(todo)[changed]
      out[i] <- y[changed]
      todo[i] <- FALSE
    }
  }
  out
}

# keys of everything that identifies an antimicrobial in `lkp` (by default the internal lookup table):
# code, name, synonyms, abbreviations, ATC codes, CID and LOINC codes
custom_ab_identifier_keys <- function(lkp = AMR_env$AB_lookup) {
  ab <- as.character(lkp$ab)
  as_list <- function(col) {
    if (!col %in% colnames(lkp)) {
      return(vector("list", length(ab)))
    }
    v <- lkp[[col]]
    if (is.list(v)) v else as.list(v)
  }
  parts <- list(as.list(ab), as_list("name"), as_list("synonyms"), as_list("abbreviations"), as_list("atc"), as_list("cid"), as_list("loinc"))
  ids <- lapply(seq_along(ab), function(i) {
    vals <- as.character(unlist(lapply(parts, function(p) p[[i]]), use.names = FALSE))
    vals[!is.na(vals) & vals != ""]
  })
  out <- data.frame(ab = rep(ab, lengths(ids)), key = custom_ab_synonym_key(unlist(ids, use.names = FALSE)), stringsAsFactors = FALSE)
  out <- out[!is.na(out$key), , drop = FALSE]
  out$pkey <- custom_ab_punct_key(out$key)
  out
}

# the same for the internal lookup table, kept until the lookup table changes (see reset_ab_cache())
custom_ab_identifier_keys_cached <- function() {
  if (is.null(AMR_env$custom_ab_id_keys)) {
    AMR_env$custom_ab_id_keys <- custom_ab_identifier_keys()
  }
  AMR_env$custom_ab_id_keys
}

# generalised names and synonyms of the internal lookup table, as used by the exact stage of as.ab()
custom_ab_generalised_ids <- function() {
  ab <- as.character(AMR_env$AB_lookup$ab)
  g <- AMR_env$AB_lookup$generalised_all
  if (!is.list(g)) g <- as.list(g)
  out <- data.frame(ab = rep(ab, lengths(g)), gen = as.character(unlist(g, use.names = FALSE)), stringsAsFactors = FALSE)
  out$genp <- custom_ab_generalised_punct(out$gen)
  out[!is.na(out$genp), , drop = FALSE]
}

#' @rdname add_custom_antimicrobials
#' @export
clear_custom_antimicrobials <- function() {
  n <- nrow(AMR_env$AB_lookup)
  AMR_env$AB_lookup <- cbind(AMR::antimicrobials, AB_LOOKUP)
  n2 <- nrow(AMR_env$AB_lookup)
  AMR_env$custom_ab_codes <- character(0)
  AMR_env$custom_ab_synonyms <- AMR_env$custom_ab_synonyms[0, , drop = FALSE]
  reset_ab_cache()
  message_("Cleared ", nr2char(n - n2), " custom record", ifelse(n - n2 > 1, "s", ""), " from the internal {.help [antimicrobials](AMR::antimicrobials)} data set.")
}

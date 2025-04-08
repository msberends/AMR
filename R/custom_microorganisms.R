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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

#' Add Custom Microorganisms
#'
#' With [add_custom_microorganisms()] you can add your own custom microorganisms, such the non-taxonomic outcome of laboratory analysis.
#' @param x A [data.frame] resembling the [microorganisms] data set, at least containing column "genus" (case-insensitive).
#' @details This function will fill in missing taxonomy for you, if specific taxonomic columns are missing, see *Examples*.
#'
#' **Important:** Due to how \R works, the [add_custom_microorganisms()] function has to be run in every \R session - added microorganisms are not stored between sessions and are thus lost when \R is exited.
#'
#' There are two ways to circumvent this and automate the process of adding microorganisms:
#'
#' **Method 1:** Using the package option [`AMR_custom_mo`][AMR-options], which is the preferred method. To use this method:
#'
#'    1. Create a data set in the structure of the [microorganisms] data set (containing at the very least column "genus") and save it with [saveRDS()] to a location of choice, e.g. `"~/my_custom_mo.rds"`, or any remote location.
#'
#'    2. Set the file location to the package option [`AMR_custom_mo`][AMR-options]: `options(AMR_custom_mo = "~/my_custom_mo.rds")`. This can even be a remote file location, such as an https URL. Since options are not saved between \R sessions, it is best to save this option to the `.Rprofile` file so that it will be loaded on start-up of \R. To do this, open the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`, add this text and save the file:
#'
#'       ```r
#'       # Add custom microorganism codes:
#'       options(AMR_custom_mo = "~/my_custom_mo.rds")
#'       ```
#'
#'       Upon package load, this file will be loaded and run through the [add_custom_microorganisms()] function.
#'
#' **Method 2:** Loading the microorganism directly from your `.Rprofile` file. Note that the definitions will be stored in a user-specific \R file, which is a suboptimal workflow. To use this method:
#'
#'    1. Edit the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`.
#'
#'    2. Add a text like below and save the file:
#'
#'       ```r
#'        # Add custom antibiotic drug codes:
#'        AMR::add_custom_microorganisms(
#'          data.frame(genus = "Enterobacter",
#'                     species = "asburiae/cloacae")
#'        )
#'       ```
#'
#' Use [clear_custom_microorganisms()] to clear the previously added microorganisms.
#' @seealso [add_custom_antimicrobials()] to add custom antimicrobials.
#' @rdname add_custom_microorganisms
#' @export
#' @examples
#' \donttest{
#' # a combination of species is not formal taxonomy, so
#' # this will result in "Enterobacter cloacae cloacae",
#' # since it resembles the input best:
#' mo_name("Enterobacter asburiae/cloacae")
#'
#' # now add a custom entry - it will be considered by as.mo() and
#' # all mo_*() functions
#' add_custom_microorganisms(
#'   data.frame(
#'     genus = "Enterobacter",
#'     species = "asburiae/cloacae"
#'   )
#' )
#'
#' # E. asburiae/cloacae is now a new microorganism:
#' mo_name("Enterobacter asburiae/cloacae")
#'
#' # its code:
#' as.mo("Enterobacter asburiae/cloacae")
#'
#' # all internal algorithms will work as well:
#' mo_name("Ent asburia cloacae")
#'
#' # and even the taxonomy was added based on the genus!
#' mo_family("E. asburiae/cloacae")
#' mo_gramstain("Enterobacter asburiae/cloacae")
#'
#' mo_info("Enterobacter asburiae/cloacae")
#'
#'
#' # the function tries to be forgiving:
#' add_custom_microorganisms(
#'   data.frame(
#'     GENUS = "BACTEROIDES / PARABACTEROIDES SLASHLINE",
#'     SPECIES = "SPECIES"
#'   )
#' )
#' mo_name("BACTEROIDES / PARABACTEROIDES")
#' mo_rank("BACTEROIDES / PARABACTEROIDES")
#'
#' # taxonomy still works, even though a slashline genus was given as input:
#' mo_family("Bacteroides/Parabacteroides")
#'
#'
#' # for groups and complexes, set them as species or subspecies:
#' add_custom_microorganisms(
#'   data.frame(
#'     genus = "Citrobacter",
#'     species = c("freundii", "braakii complex"),
#'     subspecies = c("complex", "")
#'   )
#' )
#' mo_name(c("C. freundii complex", "C. braakii complex"))
#' mo_species(c("C. freundii complex", "C. braakii complex"))
#' mo_gramstain(c("C. freundii complex", "C. braakii complex"))
#' }
add_custom_microorganisms <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  stop_ifnot("genus" %in% tolower(colnames(x)), paste0("`x` must contain column 'genus'."))

  add_MO_lookup_to_AMR_env()

  # remove any extra class/type, such as grouped tbl, or data.table:
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  colnames(x) <- tolower(colnames(x))
  # rename 'name' to 'fullname' if it's in the data set
  if ("name" %in% colnames(x) && !"fullname" %in% colnames(x)) {
    colnames(x)[colnames(x) == "name"] <- "fullname"
  }
  # keep only columns available in the microorganisms data set
  x <- x[, colnames(AMR_env$MO_lookup)[colnames(AMR_env$MO_lookup) %in% colnames(x)], drop = FALSE]

  # clean the input ----
  for (col in c("genus", "species", "subspecies")) {
    if (!col %in% colnames(x)) {
      x[, col] <- ""
    }
    if (is.factor(x[, col, drop = TRUE])) {
      x[, col] <- as.character(x[, col, drop = TRUE])
    }
    col_ <- x[, col, drop = TRUE]
    col_ <- tolower(col_)
    col_ <- gsub("slashline", "", col_, fixed = TRUE)
    col_ <- trimws2(col_)
    col_[col_ %like% "(sub)?species"] <- ""
    col_ <- gsub(" *([/-]) *", "\\1", col_, perl = TRUE)
    # groups are in our taxonomic table with a capital G
    col_ <- gsub(" group( |$)", " Group\\1", col_, perl = TRUE)

    col_[is.na(col_)] <- ""
    if (col == "genus") {
      substr(col_, 1, 1) <- toupper(substr(col_, 1, 1))
      col_ <- gsub("/([a-z])", "/\\U\\1", col_, perl = TRUE)
      stop_if(any(col_ == ""), "the 'genus' column cannot be empty")
      stop_if(any(col_ %like% " "), "the 'genus' column must not contain spaces")
    }
    x[, col] <- col_
  }
  # if subspecies is a group or complex, add it to the species and empty the subspecies
  x$species[which(x$subspecies %in% c("group", "Group", "complex"))] <- paste(
    x$species[which(x$subspecies %in% c("group", "Group", "complex"))],
    x$subspecies[which(x$subspecies %in% c("group", "Group", "complex"))]
  )
  x$subspecies[which(x$subspecies %in% c("group", "Group", "complex"))] <- ""

  if ("rank" %in% colnames(x)) {
    stop_ifnot(
      all(x$rank %in% AMR_env$MO_lookup$rank),
      "the 'rank' column can only contain these values: ", vector_or(AMR_env$MO_lookup$rank)
    )
  } else {
    x$rank <- ifelse(x$subspecies != "", "subspecies",
      ifelse(x$species != "", "species",
        ifelse(x$genus != "", "genus",
          stop("in add_custom_microorganisms(): only microorganisms up to the genus level can be added",
            call. = FALSE
          )
        )
      )
    )
  }
  x$source <- "Added by user"
  if (!"fullname" %in% colnames(x)) {
    x$fullname <- trimws2(paste(x$genus, x$species, x$subspecies))
  }
  if (!"kingdom" %in% colnames(x)) x$kingdom <- ""
  if (!"phylum" %in% colnames(x)) x$phylum <- ""
  if (!"class" %in% colnames(x)) x$class <- ""
  if (!"order" %in% colnames(x)) x$order <- ""
  if (!"family" %in% colnames(x)) x$family <- ""
  x$kingdom[is.na(x$kingdom)] <- ""
  x$phylum[is.na(x$phylum)] <- ""
  x$class[is.na(x$class)] <- ""
  x$order[is.na(x$order)] <- ""
  x$family[is.na(x$family)] <- ""

  for (col in colnames(x)) {
    if (is.factor(x[, col, drop = TRUE])) {
      x[, col] <- as.character(x[, col, drop = TRUE])
    }
    if (is.list(AMR_env$MO_lookup[, col, drop = TRUE])) {
      x[, col] <- as.list(x[, col, drop = TRUE])
    }
  }

  # fill in taxonomy based on genus
  genus_to_check <- gsub("^(.*)[^a-zA-Z].*", "\\1", x$genus, perl = TRUE)
  x$kingdom[which(x$kingdom == "" & genus_to_check != "")] <- AMR_env$MO_lookup$kingdom[match(genus_to_check[which(x$kingdom == "" & genus_to_check != "")], AMR_env$MO_lookup$genus)]
  x$phylum[which(x$phylum == "" & genus_to_check != "")] <- AMR_env$MO_lookup$phylum[match(genus_to_check[which(x$phylum == "" & genus_to_check != "")], AMR_env$MO_lookup$genus)]
  x$class[which(x$class == "" & genus_to_check != "")] <- AMR_env$MO_lookup$class[match(genus_to_check[which(x$class == "" & genus_to_check != "")], AMR_env$MO_lookup$genus)]
  x$order[which(x$order == "" & genus_to_check != "")] <- AMR_env$MO_lookup$order[match(genus_to_check[which(x$order == "" & genus_to_check != "")], AMR_env$MO_lookup$genus)]
  x$family[which(x$family == "" & genus_to_check != "")] <- AMR_env$MO_lookup$family[match(genus_to_check[which(x$family == "" & genus_to_check != "")], AMR_env$MO_lookup$genus)]

  # fill in other columns that are used in internal algorithms
  x$prevalence <- NA_real_
  x$prevalence[which(genus_to_check != "")] <- AMR_env$MO_lookup$prevalence[match(genus_to_check[which(genus_to_check != "")], AMR_env$MO_lookup$genus)]
  x$prevalence[is.na(x$prevalence)] <- 1.25
  x$status <- "accepted"
  x$ref <- paste("Self-added,", format(Sys.Date(), "%Y"))
  x$kingdom_index <- AMR_env$MO_lookup$kingdom_index[match(genus_to_check, AMR_env$MO_lookup$genus)]
  # complete missing kingdom index, so mo_matching_score() will not return NA
  x$kingdom_index[is.na(x$kingdom_index)] <- 1
  x$fullname_lower <- tolower(x$fullname)
  x$full_first <- substr(x$fullname_lower, 1, 1)
  x$species_first <- tolower(substr(x$species, 1, 1))
  x$subspecies_first <- tolower(substr(x$subspecies, 1, 1))

  if (!"mo" %in% colnames(x)) {
    # create the mo code
    x$mo <- NA_character_
  }
  x$mo <- trimws2(as.character(x$mo))
  x$mo[x$mo == ""] <- NA_character_
  current <- sum(AMR_env$MO_lookup$source == "Added by user", na.rm = TRUE)
  x$mo[is.na(x$mo)] <- paste0(
    "CUSTOM",
    seq.int(from = current + 1, to = current + nrow(x), by = 1),
    "_",
    trimws(
      paste(abbreviate_mo(x$genus, 5),
        abbreviate_mo(x$species, 4, hyphen_as_space = TRUE),
        abbreviate_mo(x$subspecies, 4, hyphen_as_space = TRUE),
        sep = "_"
      ),
      whitespace = "_"
    )
  )
  stop_if(anyDuplicated(c(as.character(AMR_env$MO_lookup$mo), x$mo)), "MO codes must be unique and not match existing MO codes of the AMR package")

  # add to package ----
  AMR_env$custom_mo_codes <- c(AMR_env$custom_mo_codes, x$mo)
  class(AMR_env$MO_lookup$mo) <- "character"

  new_df <- AMR_env$MO_lookup[0, , drop = FALSE][seq_len(NROW(x)), , drop = FALSE]
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

  # clear previous coercions
  suppressMessages(mo_reset_session())

  AMR_env$MO_lookup <- unique(rbind_AMR(AMR_env$MO_lookup, new_df))
  class(AMR_env$MO_lookup$mo) <- c("mo", "character")
  if (nrow(x) <= 3) {
    message_("Added ", vector_and(italicise(x$fullname), quotes = FALSE), " to the internal `microorganisms` data set.")
  } else {
    message_("Added ", nr2char(nrow(x)), " records to the internal `microorganisms` data set.")
  }
}

#' @rdname add_custom_microorganisms
#' @export
clear_custom_microorganisms <- function() {
  n <- nrow(AMR_env$MO_lookup)

  # reset
  AMR_env$MO_lookup <- NULL
  add_MO_lookup_to_AMR_env()

  # clear previous coercions
  suppressMessages(mo_reset_session())

  n2 <- nrow(AMR_env$MO_lookup)
  AMR_env$custom_mo_codes <- character(0)
  AMR_env$mo_previously_coerced <- AMR_env$mo_previously_coerced[which(AMR_env$mo_previously_coerced$mo %in% AMR_env$MO_lookup$mo), , drop = FALSE]
  AMR_env$mo_uncertainties <- AMR_env$mo_uncertainties[0, , drop = FALSE]
  message_("Cleared ", nr2char(n - n2), " custom record", ifelse(n - n2 > 1, "s", ""), " from the internal `microorganisms` data set.")
}

abbreviate_mo <- function(x, minlength = 5, prefix = "", hyphen_as_space = FALSE, ...) {
  if (hyphen_as_space == TRUE) {
    x <- gsub("-", " ", x, fixed = TRUE)
  }
  # keep a starting Latin ae
  suppressWarnings(
    gsub(
      "(\u00C6|\u00E6)+",
      "AE",
      toupper(
        paste0(
          prefix,
          abbreviate(
            gsub("^ae",
              "\u00E6\u00E6",
              x,
              ignore.case = TRUE
            ),
            minlength = minlength,
            use.classes = TRUE,
            method = "both.sides",
            ...
          )
        )
      )
    )
  )
}

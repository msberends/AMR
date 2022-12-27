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

#' Add Custom Microorganisms to This Package
#'
#' With [add_custom_microorganisms()] you can add your own custom antimicrobial drug codes to the `AMR` package.
#' @param x a [data.frame] resembling the [microorganisms] data set, at least containing columns "mo", "genus" and "species"
#' @details This function will fill in missing taxonomy for you, if specific taxonomic columns are missing, see *Examples*.
#' 
#' **Important:** Due to how \R works, the [add_custom_microorganisms()] function has to be run in every \R session - added microorganisms are not stored between sessions and are thus lost when \R is exited. 
#' 
#' There are two ways to automate this process:
#' 
#' **Method 1:** Save the microorganisms to a local or remote file (can even be the internet). To use this method:
#' 
#'    1. Create a data set in the structure of the [microorganisms] data set (containing at the very least columns "ab" and "name") and save it with [saveRDS()] to a location of choice, e.g. `"~/my_custom_mo.rds"`, or any remote location.
#'    
#'    2. Set the file location to the `AMR_custom_mo` \R option: `options(AMR_custom_mo = "~/my_custom_mo.rds")`. This can even be a remote file location, such as an https URL. Since options are not saved between \R sessions, it is best to save this option to the `.Rprofile` file so that it will loaded on start-up of \R. To do this, open the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`, add this text and save the file:
#'
#'       ```r
#'       # Add custom microorganism codes:
#'       options(AMR_custom_mo = "~/my_custom_mo.rds")
#'       ```
#'       
#'       Upon package load, this file will be loaded and run through the [add_custom_microorganisms()] function.
#' 
#' **Method 2:** Save the microorganism directly to your `.Rprofile` file. An important downside is that this requires to load the `AMR` package at every start-up. To use this method:
#' 
#'    1. Edit the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`.
#'
#'    2. Add a text like below and save the file:
#'
#'       ```r
#'        # Add custom antibiotic drug codes:
#'        library(AMR)
#'        add_custom_microorganisms(
#'          data.frame(mo = "ENT_ASB_CLO",
#'                     genus = "Enterobacter",
#'                     species = "asburiae/cloacae")
#'        )
#'       ```
#'
#' Use [clear_custom_microorganisms()] to clear the previously added antimicrobials.
#' @seealso [add_custom_antimicrobials()] to add custom antimicrobials to this package.
#' @rdname add_custom_microorganisms
#' @export
#' @examples
#' \donttest{
#'
#' # a combination of species is not formal taxonomy, so
#' # this will result in only "Enterobacter asburiae":
#' mo_name("Enterobacter asburiae/cloacae")
#'
#' # now add a custom entry - it will be considered by as.mo() and
#' # all mo_*() functions
#' add_custom_microorganisms(
#'   data.frame(mo = "ENT_ASB_CLO",
#'              genus = "Enterobacter",
#'              species = "asburiae/cloacae"
#'   )
#' )
#'
#' # "ENT_ASB_CLO" is now a new microorganism:
#' mo_name("Enterobacter asburiae/cloacae")
#' as.mo("ent_asb_clo")
#' mo_name("ent_asb_clo")
#' # all internal algorithms will work as well:
#' mo_name("Ent asburia cloacae")
#' 
#' # and even the taxonomy was added based on the genus!
#' mo_family("ent_asb_clo")
#' mo_gramstain("Enterobacter asburiae/cloacae")
#'
#' mo_info("ent_asb_clo")
#' }
add_custom_microorganisms <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  required_cols <- c("mo", "genus", "species")
  stop_ifnot(
    all(required_cols %in% colnames(x)),
    paste0("`x` must contain columns ", vector_and(required_cols, sort = FALSE), ".")
  )
  stop_if(
    any(x$mo %in% AMR_env$MO_lookup$mo),
    "Microorganism code(s) ", vector_and(x$mo[x$mo %in% AMR_env$MO_lookup$mo]), " already exist in the internal `microorganisms` data set."
  )
  # remove any extra class/type, such as grouped tbl, or data.table:
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  # rename 'name' to 'fullname' if it's in the data set
  if ("name" %in% colnames(x) && !"fullname" %in% colnames(x)) {
    colnames(x)[colnames(x) == "name"] <- "fullname"
  }
  # keep only columns available in the microorganisms data set
  x <- x[, colnames(AMR_env$MO_lookup)[colnames(AMR_env$MO_lookup) %in% colnames(x)], drop = FALSE]
  
  # clean the input ----
  if (!"subspecies" %in% colnames(x)) {
    x$subspecies <- NA_character_
  }
  x$genus <- trimws2(x$genus)
  x$species <- trimws2(x$species)
  x$subspecies <- trimws2(x$subspecies)
  x$genus[x$genus == ""] <- NA_character_
  x$species[x$species == ""] <- NA_character_
  x$subspecies[x$subspecies == ""] <- NA_character_
  stop_if(any(x$genus[!is.na(x$genus)] %like% " "),
          "the 'genus' column must not contain spaces")
  stop_if(any(x$species[!is.na(x$species)] %like% " "),
          "the 'species' column must not contain spaces")
  stop_if(any(x$subspecies[!is.na(x$subspecies)] %like% " "),
          "the 'subspecies' column must not contain spaces")
  
  if ("rank" %in% colnames(x)) {
    stop_ifnot(all(x$rank %in% AMR_env$MO_lookup$rank),
               "the 'rank' column can only contain these values: ", vector_or(AMR_env$MO_lookup$rank))
  } else {
    x$rank <- ifelse(!is.na(x$subspecies), "subspecies",
                     ifelse(!is.na(x$species), "species",
                            ifelse(!is.na(x$genus), "genus",
                                   stop("in add_custom_microorganisms(): the 'genus' column cannot be empty",
                                        call. = FALSE))))
  }
  if (!"fullname" %in% colnames(x)) {
    x$fullname <- paste(x$genus, x$species, x$subspecies)
    x$fullname <- gsub(" NA", "", x$fullname)
  }
  if (!"kingdom" %in% colnames(x)) x$kingdom <- NA_character_
  if (!"phylum" %in% colnames(x)) x$phylum <- NA_character_
  if (!"class" %in% colnames(x)) x$class <- NA_character_
  if (!"order" %in% colnames(x)) x$order <- NA_character_
  if (!"family" %in% colnames(x)) x$family <- NA_character_

  for (col in colnames(x)) {
    if (is.list(AMR_env$MO_lookup[, col, drop = TRUE])) {
      x[, col] <- as.list(x[, col, drop = TRUE])
    }
  }
  
  # fill in other columns
  x$status <- "accepted"
  x$prevalence <- 1
  x$kingdom <- AMR_env$MO_lookup$kingdom[match(x$genus, AMR_env$MO_lookup$genus)]
  x$phylum <- AMR_env$MO_lookup$phylum[match(x$genus, AMR_env$MO_lookup$genus)]
  x$class <- AMR_env$MO_lookup$class[match(x$genus, AMR_env$MO_lookup$genus)]
  x$order <- AMR_env$MO_lookup$order[match(x$genus, AMR_env$MO_lookup$genus)]
  x$family <- AMR_env$MO_lookup$family[match(x$genus, AMR_env$MO_lookup$genus)]
  
  x$kingdom_index <- AMR_env$MO_lookup$kingdom_index[match(x$genus, AMR_env$MO_lookup$genus)]
  x$fullname_lower <- tolower(x$fullname)
  x$full_first <- substr(x$fullname_lower, 1, 1)
  x$species_first <- tolower(substr(x$species, 1, 1))
  x$subspecies_first <- tolower(substr(x$subspecies, 1, 1))
  
  # add to pacakge ----
  
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
  AMR_env$MO_lookup <- unique(rbind(AMR_env$MO_lookup, new_df))
  AMR_env$mo_previously_coerced <- AMR_env$mo_previously_coerced[which(!AMR_env$mo_previously_coerced$mo %in% new_df$mo), , drop = FALSE]
  class(AMR_env$MO_lookup$mo) <- c("mo", "character")
  message_("Added ", nr2char(nrow(x)), " record", ifelse(nrow(x) > 1, "s", ""), " to the internal `microorganisms` data set.")
}

#' @rdname add_custom_microorganisms
#' @export
clear_custom_microorganisms <- function() {
  n <- nrow(AMR_env$MO_lookup)
  AMR_env$MO_lookup <- create_MO_lookup()
  n2 <- nrow(AMR_env$MO_lookup)
  AMR_env$custom_mo_codes <- character(0)
  AMR_env$mo_previously_coerced <- AMR_env$mo_previously_coerced[which(AMR_env$mo_previously_coerced$mo %in% AMR_env$MO_lookup$mo), , drop = FALSE]
  AMR_env$mo_uncertainties <- AMR_env$mo_uncertainties[0, , drop = FALSE]
  message_("Cleared ", nr2char(n - n2), " custom record", ifelse(n - n2 > 1, "s", ""), " from the internal `microorganisms` data set.")
}

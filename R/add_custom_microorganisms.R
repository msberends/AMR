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
#' With [add_custom_microorganisms()] you can add your own custom microorganisms to the `AMR` package, such the non-taxonomic outcome of laboratory analysis.
#' @param x a [data.frame] resembling the [microorganisms] data set, at least containing columns "genus" and "species"
#' @details This function will fill in missing taxonomy for you, if specific taxonomic columns are missing, see *Examples*.
#' 
#' **Important:** Due to how \R works, the [add_custom_microorganisms()] function has to be run in every \R session - added microorganisms are not stored between sessions and are thus lost when \R is exited. 
#' 
#' There are two ways to automate this process:
#' 
#' **Method 1:** Save the microorganisms to a local or remote file (can even be the internet). To use this method:
#' 
#'    1. Create a data set in the structure of the [microorganisms] data set (containing at the very least columns "genus" and "species") and save it with [saveRDS()] to a location of choice, e.g. `"~/my_custom_mo.rds"`, or any remote location.
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
#'          data.frame(genus = "Enterobacter",
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
#'   data.frame(genus = "Enterobacter",
#'              species = "asburiae/cloacae"
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
#' }
add_custom_microorganisms <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  stop_ifnot(
    all(c("genus", "species") %in% colnames(x)),
    paste0("`x` must contain columns ", vector_and(c("genus", "species"), sort = FALSE), ".")
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
  x$genus[is.na(x$genus)] <- ""
  x$species[is.na(x$species)] <- ""
  x$subspecies[is.na(x$subspecies)] <- ""
  stop_if(any(x$genus %like% " "),
          "the 'genus' column must not contain spaces")
  stop_if(any(x$species %like% " "),
          "the 'species' column must not contain spaces")
  stop_if(any(x$subspecies %like% " "),
          "the 'subspecies' column must not contain spaces")
  
  if ("rank" %in% colnames(x)) {
    stop_ifnot(all(x$rank %in% AMR_env$MO_lookup$rank),
               "the 'rank' column can only contain these values: ", vector_or(AMR_env$MO_lookup$rank))
  } else {
    x$rank <- ifelse(x$subspecies != "", "subspecies",
                     ifelse(x$species != "", "species",
                            ifelse(x$genus != "", "genus",
                                   stop("in add_custom_microorganisms(): the 'genus' column cannot be empty",
                                        call. = FALSE))))
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
    if (is.list(AMR_env$MO_lookup[, col, drop = TRUE])) {
      x[, col] <- as.list(x[, col, drop = TRUE])
    }
  }
  
  # fill in other columns
  x$status <- "accepted"
  x$prevalence <- 1
  
  x$kingdom[which(x$kingdom == "" & x$genus != "")] <- AMR_env$MO_lookup$kingdom[match(x$genus[which(x$kingdom == "" & x$genus != "")], AMR_env$MO_lookup$genus)]
  x$phylum[which(x$phylum == "" & x$genus != "")] <- AMR_env$MO_lookup$phylum[match(x$genus[which(x$phylum == "" & x$genus != "")], AMR_env$MO_lookup$genus)]
  x$class[which(x$class == "" & x$genus != "")] <- AMR_env$MO_lookup$class[match(x$genus[which(x$class == "" & x$genus != "")], AMR_env$MO_lookup$genus)]
  x$order[which(x$order == "" & x$genus != "")] <- AMR_env$MO_lookup$order[match(x$genus[which(x$order == "" & x$genus != "")], AMR_env$MO_lookup$genus)]
  x$family[which(x$family == "" & x$genus != "")] <- AMR_env$MO_lookup$family[match(x$genus[which(x$family == "" & x$genus != "")], AMR_env$MO_lookup$genus)]
  
  x$kingdom_index <- AMR_env$MO_lookup$kingdom_index[match(x$genus, AMR_env$MO_lookup$genus)]
  x$fullname_lower <- tolower(x$fullname)
  x$full_first <- substr(x$fullname_lower, 1, 1)
  x$species_first <- tolower(substr(x$species, 1, 1))
  x$subspecies_first <- tolower(substr(x$subspecies, 1, 1))
  
  if (!"mo" %in% colnames(x)) {
    # create the mo code
    x$mo <- NA_character_
  }
  x$mo <- trimws2(x$mo)
  x$mo[x$mo == ""] <- NA_character_
  x$mo[is.na(x$mo)] <- paste0("CUSTOM",
                              seq(from = sum(AMR_env$MO_lookup$source == "Added by user", na.rm = TRUE) + 1, to = nrow(x), by = 1),
                              "_",
                              toupper(unname(abbreviate(gsub(" +", " _ ",
                                                             gsub("[^A-Za-z0-9-]", " ",
                                                                  trimws2(paste(x$genus, x$species, x$subspecies)))),
                                                        minlength = 10))))
  
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
  
  AMR_env$MO_lookup <- unique(rbind(AMR_env$MO_lookup, new_df))
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

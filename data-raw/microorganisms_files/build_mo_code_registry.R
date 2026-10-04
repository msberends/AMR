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

# This script (re)creates data-raw/microorganisms_files/mo_code_registry.csv: the registry of all MO codes
# that were ever part of a released version of the AMR package, since v2.0.0. See the README.md in this
# folder for the rules.
#
# Run it from the root of the repository, after every new release (once its git tag exists):
#   source("data-raw/microorganisms_files/build_mo_code_registry.R")
# It always rebuilds the whole registry from the git tags, so it is reproducible and can be run at any time.
# Development versions (the main branch) are never part of the registry.

library(dplyr)
devtools::load_all(".", quiet = TRUE) # for %like_case%

registry_file <- "data-raw/microorganisms_files/mo_code_registry.csv"

# Decision by Matthijs S. Berends, 3 October 2026: MO codes from releases before v2.0.0 (12 March 2023)
# are not supported.
first_supported_release <- "2.0.0"

# all release tags (vX.Y.Z, so no development versions) since the first supported release
tags <- system2("git", "tag", stdout = TRUE)
tags <- tags[tags %like_case% "^v[0-9]+[.][0-9]+[.][0-9]+$"]
tags <- tags[numeric_version(sub("^v", "", tags)) >= numeric_version(first_supported_release)]
tags <- tags[order(numeric_version(sub("^v", "", tags)))]
message("Releases: ", toString(tags))

read_release <- function(tag) {
  tmp <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp))
  status <- system2("git", c("show", paste0(tag, ":data/microorganisms.rda")), stdout = tmp)
  if (!identical(status, 0L)) {
    stop("Could not read data/microorganisms.rda of ", tag, call. = FALSE)
  }
  env <- new.env()
  load(tmp, envir = env)
  mo <- as.data.frame(env$microorganisms, stringsAsFactors = FALSE)
  # until v3.0.1, `kingdom` contained what is now called `domain`
  if (!"domain" %in% colnames(mo)) {
    mo$domain <- mo$kingdom
    mo$domain[mo$domain == "(unknown kingdom)"] <- "(unknown domain)"
  }
  data.frame(
    release = tag,
    mo = as.character(mo$mo),
    fullname = mo$fullname,
    rank = mo$rank,
    domain = mo$domain,
    genus = mo$genus,
    species = mo$species,
    subspecies = mo$subspecies,
    stringsAsFactors = FALSE
  )
}

all_releases <- bind_rows(lapply(tags, read_release)) %>%
  mutate(release_order = match(release, tags))

# one row per MO code; the name, rank and taxonomy are taken from the most recent release that contained the code
registry <- all_releases %>%
  group_by(mo) %>%
  arrange(desc(release_order), .by_group = TRUE) %>%
  summarise(
    # other names the code had in older releases (e.g. a corrected subspecies name), for transparency
    # (must come first, since `fullname` is overwritten below)
    previous_names = paste(setdiff(unique(fullname), first(fullname)), collapse = "; "),
    fullname = first(fullname),
    rank = first(rank),
    domain = first(domain),
    genus = first(genus),
    species = first(species),
    subspecies = first(subspecies),
    first_release = tags[min(release_order)],
    last_release = tags[max(release_order)],
    .groups = "drop"
  ) %>%
  relocate(previous_names, .after = everything()) %>%
  arrange(mo)

message(
  nrow(registry), " MO codes, of which ", sum(registry$previous_names != ""),
  " had another name in an older release (check these when the registry changes)"
)

write.csv(registry, registry_file, row.names = FALSE, na = "")
message("Saved to ", registry_file)

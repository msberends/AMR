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

# Build the function and data set reference of one released AMR version as Markdown only, for use by LLMs, at https://amr-for-r.org/release/<version>/.
# Used by .github/workflows/website-release.yaml, after installing AMR from the release tag:
#   Rscript build_release_reference.R <package source of the tag> <checkout of the gh-pages branch> <version>
# Writes <gh-pages>/release/<version>/reference/*.md and <gh-pages>/release/<version>/llms.txt, and regenerates the version index <gh-pages>/release/llms.txt.
# The HTML pages that pkgdown renders are only an intermediate step and are not published, so the main website remains the only human-readable reference.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript build_release_reference.R <package source> <gh-pages checkout> <version>")
}
src <- normalizePath(args[1], mustWork = TRUE)
site <- normalizePath(args[2], mustWork = TRUE)
version <- args[3]
if (!grepl("^[0-9]+[.][0-9]+[.][0-9]+$", version)) {
  stop("Version must be a released version such as 3.1.0, not '", version, "'")
}

desc <- read.dcf(file.path(src, "DESCRIPTION"))
if (desc[1, "Version"] != version) {
  stop("The package source has version ", desc[1, "Version"], ", not ", version)
}
release_date <- if ("Date" %in% colnames(desc)) desc[1, "Date"] else NA_character_

base_url <- paste0("https://amr-for-r.org/release/", version)
html_dir <- file.path(tempfile("pkgdown_"), version)
release_root <- file.path(site, "release")
out_dir <- file.path(release_root, version)

# render the reference as HTML in a temporary folder, without running examples
pkg <- pkgdown::as_pkgdown(src, override = list(destination = html_dir, url = base_url))
pkgdown::build_reference(pkg, examples = FALSE, lazy = FALSE, devel = FALSE)

# convert every HTML page to Markdown, with absolute links to the Markdown pages of this same version
unlink(out_dir, recursive = TRUE)
dir.create(file.path(out_dir, "reference"), recursive = TRUE)
html_files <- list.files(file.path(html_dir, "reference"), pattern = "[.]html$")
for (html_file in html_files) {
  pkgdown:::convert_md(
    src_path = file.path(html_dir, "reference", html_file),
    dst_path = file.path(out_dir, "reference", sub("[.]html$", ".md", html_file)),
    url = paste0(base_url, "/reference/", html_file)
  )
}
md_files <- list.files(file.path(out_dir, "reference"), pattern = "[.]md$")
if (length(md_files) != length(html_files) || !"index.md" %in% md_files) {
  stop("Converting HTML to Markdown failed: ", length(html_files), " HTML pages, ", length(md_files), " Markdown pages")
}
unlink(dirname(html_dir), recursive = TRUE)

# llms.txt of this version: a header that states the version, followed by the reference index
index_md <- readLines(file.path(out_dir, "reference", "index.md"), encoding = "UTF-8")
index_md <- index_md[-1] # drop the '# Package index' title, the header below replaces it
index_md <- index_md[cumsum(index_md != "") > 0] # and the blank lines that followed it
writeLines(
  c(
    paste0("# AMR for R ", version, ": reference of all functions and data sets"),
    "",
    paste0(
      "> This is the complete reference documentation of version ", version, " of the AMR package for R",
      if (!is.na(release_date)) paste0(", released on ", release_date) else "",
      ". It applies to AMR ", version, " only: functions, arguments, defaults, data sets and their columns can differ in other versions."
    ),
    "",
    "How to use this reference:",
    "",
    paste0("- Every function and data set has its own page, linked below. The pages are Markdown and live at `", base_url, "/reference/<page>.md`."),
    "- Each page contains the usage (all arguments with their defaults), the description of every argument, details, the return value and examples.",
    "- Several related functions share one page; the list below shows which page documents which function (e.g., `mo_gramstain()` is documented on `mo_property.md`).",
    "- A function or argument that is not in this reference does not exist in this version.",
    "- The documentation of other versions is listed at <https://amr-for-r.org/release/llms.txt>.",
    "",
    index_md
  ),
  file.path(out_dir, "llms.txt"),
  useBytes = TRUE
)

# regenerate the version index from all versions present in the gh-pages checkout
versions <- list.files(release_root, pattern = "^[0-9]+[.][0-9]+[.][0-9]+$")
versions <- versions[file.exists(file.path(release_root, versions, "llms.txt"))]
versions <- versions[order(numeric_version(versions), decreasing = TRUE)]
version_dates <- vapply(
  versions,
  function(v) {
    first_lines <- readLines(file.path(release_root, v, "llms.txt"), n = 3, encoding = "UTF-8")
    date <- regmatches(first_lines, regexpr("released on [0-9]{4}-[0-9]{2}-[0-9]{2}", first_lines))
    if (length(date) == 0) "" else paste0(", ", date[1])
  },
  character(1)
)
writeLines(
  c(
    "# AMR for R: reference documentation per package version",
    "",
    "> The AMR package for R changes between versions: functions, arguments, defaults, data sets and clinical breakpoints are added or changed. Answers about AMR must therefore be based on the documentation of the version that the user has installed.",
    "",
    "How to find the right documentation:",
    "",
    "1. Determine the installed version, which is the output of `packageVersion(\"AMR\")` in R.",
    "2. If the version has three numbers (e.g., `3.1.0`), it is a release from CRAN. Use the reference of exactly that version from the list below.",
    "3. If the version has four numbers (e.g., `3.1.0.9022`), it is a development version from GitHub. Use the reference of the latest development version at <https://amr-for-r.org/llms.txt>.",
    "4. If the installed version is not listed below, use the closest older version and state this.",
    "",
    "## Released versions",
    "",
    paste0("- [AMR ", versions, "](https://amr-for-r.org/release/", versions, "/llms.txt): reference of all functions and data sets", version_dates)
  ),
  file.path(release_root, "llms.txt"),
  useBytes = TRUE
)

message("Built the reference of AMR ", version, ": ", length(md_files), " Markdown pages in ", out_dir)

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

# Run R CMD check on all CRAN packages that depend on AMR, against the INSTALLED development version of AMR.
# Used by .github/workflows/check-revdeps.yaml; run from the root of the repository, after installing AMR from it:
#   Rscript .github/scripts/check_revdeps.R
# Fails if any reverse dependency returns an ERROR or WARNING.

# only CRAN, also when other repositories are set, since this is about CRAN releases
repos <- c(CRAN = "https://cloud.r-project.org")
db <- utils::available.packages(repos = repos)
revdeps <- sort(tools::package_dependencies("AMR", db = db, reverse = TRUE, which = "all")[["AMR"]])
if (length(revdeps) == 0) {
  message("No reverse dependencies on CRAN")
  quit(status = 0)
}
message("Reverse dependencies: ", paste(revdeps, collapse = ", "))

# all dependencies of the reverse dependencies, but never AMR itself, which is the development version installed from this repository
deps <- unique(unlist(tools::package_dependencies(revdeps, db = db, which = "all"), use.names = FALSE))
base_pkgs <- rownames(utils::installed.packages(priority = "base"))
deps <- setdiff(deps, c("AMR", base_pkgs, rownames(utils::installed.packages())))
if (length(deps) > 0) {
  # pak also installs the system requirements of these packages
  tryCatch(pak::pkg_install(deps, upgrade = FALSE), error = function(e) {
    message("pak failed (", conditionMessage(e), "), falling back to install.packages()")
    utils::install.packages(deps)
  })
}
amr_version <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
if (as.character(utils::packageVersion("AMR")) != amr_version) {
  stop("installed AMR is v", utils::packageVersion("AMR"), ", expected the development version v", amr_version)
}
message("Checking against AMR v", amr_version)

check_dir <- file.path(tempdir(), "revdeps")
dir.create(check_dir, showWarnings = FALSE)
tarballs <- utils::download.packages(revdeps, destdir = check_dir, repos = repos, type = "source")[, 2]

results <- vapply(tarballs, FUN.VALUE = character(1), function(tarball) {
  pkg <- sub("_.*", "", basename(tarball))
  message("\n=== Checking ", pkg, " ===")
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "check", "--no-manual", "--as-cran", "-o", shQuote(check_dir), shQuote(tarball)),
    stdout = TRUE,
    stderr = TRUE,
    # skip the incoming checks for CRAN submission, which concern the maintainers of the reverse dependencies
    env = c("_R_CHECK_FORCE_SUGGESTS_=false", "_R_CHECK_CRAN_INCOMING_=false")
  ))
  cat(out, sep = "\n")
  status <- grep("^Status:", out, value = TRUE)
  if (length(status) == 0) "ERROR (check did not finish)" else trimws(sub("^Status:", "", status[length(status)]))
}, USE.NAMES = FALSE)
names(results) <- revdeps

message("\n=== Summary against AMR v", amr_version, " ===")
message(paste0(format(names(results)), "  ", results, collapse = "\n"))
if (any(grepl("ERROR|WARNING", results))) {
  quit(status = 1)
}

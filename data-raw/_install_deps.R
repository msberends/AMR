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

# some old R instances have trouble installing tinytest, so we ship it too
install.packages("data-raw/tinytest_1.3.1.tar.gz", dependencies = c("Depends", "Imports", "LinkingTo"))
install.packages(getwd(), repos = NULL, type = "source")

pkg_suggests <- gsub(
  "[^a-zA-Z0-9.]+", "",
  unlist(strsplit(unlist(packageDescription("AMR",
    fields = c("Suggests", "Enhances", "LinkingTo")
  )),
  split = ", ?"
  ))
)
pkg_suggests <- unname(pkg_suggests[!is.na(pkg_suggests)])
cat("################################################\n")
cat("Packages listed in Suggests/Enhances:", paste(pkg_suggests, collapse = ", "), "\n")
cat("################################################\n")

if (.Platform$OS.type != "unix") {
  # no compiling on Windows here
  options(install.packages.compile.from.source = FALSE)
}

to_install <- pkg_suggests[!pkg_suggests %in% rownames(utils::installed.packages())]
if (length(to_install) == 0) {
  message("\nNothing to install\n")
}
for (i in seq_len(length(to_install))) {
  cat("Installing package", to_install[i], "\n")
  tryCatch(install.packages(to_install[i],
    type = "source",
    repos = "https://cran.rstudio.com/",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    quiet = FALSE
  ),
  # message = function(m) invisible(),
  warning = function(w) message(w$message),
  error = function(e) message(e$message)
  )
  if (.Platform$OS.type != "unix" && !to_install[i] %in% rownames(utils::installed.packages())) {
    tryCatch(install.packages(to_install[i],
      type = "binary",
      repos = "https://cran.rstudio.com/",
      dependencies = c("Depends", "Imports", "LinkingTo"),
      quiet = FALSE
    ),
    # message = function(m) invisible(),
    warning = function(w) message(w$message),
    error = function(e) message(e$message)
    )
  }
}

to_update <- as.data.frame(utils::old.packages(repos = "https://cran.rstudio.com/"), stringsAsFactors = FALSE)
to_update <- to_update[which(to_update$Package %in% pkg_suggests), "Package", drop = TRUE]
if (length(to_update) == 0) {
  message("\nNothing to update\n")
}
for (i in seq_len(length(to_update))) {
  cat("Updating package '", to_update[i], "' v", as.character(packageVersion(to_update[i])), "\n", sep = "")
  tryCatch(update.packages(to_update[i], repos = "https://cran.rstudio.com/", ask = FALSE),
    # message = function(m) invisible(),
    warning = function(w) message(w$message),
    error = function(e) message(e$message)
  )
  cat("Updated to '", to_update[i], "' v", as.character(packageVersion(to_update[i])), "\n", sep = "")
}

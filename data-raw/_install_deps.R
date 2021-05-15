# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
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

pkg_suggests <- AMR:::trimws(unlist(strsplit(packageDescription("AMR")$Suggests, ",(\n)?")))

to_install <- pkg_suggests[!pkg_suggests %in% rownames(utils::installed.packages())]
to_update <- as.data.frame(old.packages(), stringsAsFactors = FALSE)

for (i in seq_len(length(to_install))) {
  cat("Installing package", to_install[i], "\n")
  tryCatch(install.packages(to_install[i], repos = "https://cran.rstudio.com/", dependencies = TRUE, quiet = TRUE),
           message = function(m) invisible(),
           warning = function(w) message(w$message),
           error = function(e) message(e$message))
}

for (i in seq_len(length(to_update))) {
  cat("Updating package", to_install[i], "\n")
  tryCatch(update.packages(to_update[i], repos = "https://cran.rstudio.com/", ask = FALSE),
           message = function(m) invisible(),
           warning = function(w) message(w$message),
           error = function(e) message(e$message))
}

# saveRDS(to_update, ".github/depends.Rds", version = 2)

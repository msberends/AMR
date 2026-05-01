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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Download and Unpack an AMR Course Repository
#'
#' Downloads and unpacks a GitHub repository containing course materials, using [usethis::use_course()]. This is a convenience wrapper intended for use in educational settings, such as workshops or tutorials associated with the AMR package.
#' @param github_repo A character string specifying the GitHub repository with username and repo name, e.g. `"https://github.com/username/repo"`.
#' @param branch A character string specifying the branch to download. Defaults to `"main"`.
#' @param ... Additional arguments passed on to [usethis::use_course()].
#' @details
#' This function constructs a ZIP archive URL from the provided `github_repo` and `branch`, then delegates to [usethis::use_course()] to handle the download and extraction.
#'
#' The function is designed for interactive use in course or workshop settings and is not intended for use in non-interactive or automated pipelines.
#' @return
#' Called for its side effect. [usethis::use_course()] will prompt the user to choose a destination and open the extracted project. Returns invisibly whatever [usethis::use_course()] returns.
#' @seealso [usethis::use_course()]
#' @export
#' @examples
#' \dontrun{
#'
#' # Let this run by users, e.g., webinar participants
#' amr_course("https://github.com/my_user_name/our_AMR_course")
#' }
amr_course <- function(github_repo, branch = "main", ...) {
  if (!"usethis" %in% rownames(utils::installed.packages())) {
    if ("rlang" %in% rownames(utils::installed.packages())) {
      rlang::check_installed("usethis")
    } else {
      stop("Package usethis is not installed. Please run: install.packages(\"usethis\")", call. = FALSE)
    }
  }
  url <- paste0(github_repo, "/archive/refs/heads/", branch, ".zip")
  use_course <- import_fn("use_course", "usethis")
  message("This will download and unpack the contents of a repository.\n")
  use_course(url, ...)
}

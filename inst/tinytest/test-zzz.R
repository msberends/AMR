# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# Check if these functions still exist in their package (all are in Suggests field)
# Since GitHub Actions runs every night, we will get emailed when a dependency fails based on this unit test

# functions used by import_fn()
import_functions <- c(
  "%chin%" = "data.table",
  "anti_join" = "dplyr",
  "chmatch" = "data.table",
  "cur_column" = "dplyr",
  "full_join" = "dplyr",
  "has_internet" = "curl",
  "html_attr" = "rvest",
  "html_children" = "rvest",
  "html_node" = "rvest",
  "html_nodes" = "rvest",
  "html_table" = "rvest",
  "html_text" = "rvest",
  "inner_join" = "dplyr",
  "insertText" = "rstudioapi",
  "left_join" = "dplyr",
  "new_pillar_shaft_simple" = "pillar",
  "progress_bar" = "progress",
  "read_html" = "xml2",
  "right_join" = "dplyr",
  "semi_join" = "dplyr",
  "showQuestion" = "rstudioapi"
)

# functions that are called directly with ::
call_functions <- c(
  # cleaner
  "freq.default" = "cleaner",
  # cli
  "symbol" = "cli",
  "ansi_has_hyperlink_support" = "cli",
  # rstudioapi (RStudio)
  "isAvailable" = "rstudioapi",
  "versionInfo" = "rstudioapi",
  # readxl
  "read_excel" = "readxl",
  # ggplot2
  "aes" = "ggplot2",
  "aes_string" = "ggplot2",
  "arrow" = "ggplot2",
  "autoplot" = "ggplot2",
  "element_blank" = "ggplot2",
  "element_line" = "ggplot2",
  "element_text" = "ggplot2",
  "expand_limits" = "ggplot2",
  "facet_wrap" = "ggplot2",
  "geom_errorbar" = "ggplot2",
  "geom_path" = "ggplot2",
  "geom_point" = "ggplot2",
  "geom_ribbon" = "ggplot2",
  "geom_segment" = "ggplot2",
  "geom_text" = "ggplot2",
  "ggplot" = "ggplot2",
  "labs" = "ggplot2",
  "layer" = "ggplot2",
  "position_fill" = "ggplot2",
  "position_dodge2" = "ggplot2",
  "scale_fill_manual" = "ggplot2",
  "scale_y_continuous" = "ggplot2",
  "theme" = "ggplot2",
  "theme_minimal" = "ggplot2",
  "unit" = "ggplot2",
  "xlab" = "ggplot2",
  "ylab" = "ggplot2"
)
if (AMR:::pkg_is_available("skimr", min_version = "2.0.0")) {
  call_functions <- c(call_functions,
    # skimr
    "inline_hist" = "skimr",
    "sfl" = "skimr"
  )
}

extended_functions <- c(
  "freq" = "cleaner",
  "autoplot" = "ggplot2",
  "pillar_shaft" = "pillar",
  "get_skimmers" = "skimr",
  "type_sum" = "pillar",
  "vec_cast" = "vctrs",
  "vec_math" = "vctrs",
  "vec_ptype2" = "vctrs"
)

import_functions <- c(import_functions, call_functions, extended_functions)
for (i in seq_len(length(import_functions))) {
  fn <- names(import_functions)[i]
  pkg <- unname(import_functions[i])
  # function should exist in foreign pkg namespace
  if (AMR:::pkg_is_available(pkg,
    also_load = FALSE,
    min_version = if (pkg == "dplyr") "1.0.0" else NULL
  )) {
    expect_true(!is.null(AMR:::import_fn(name = fn, pkg = pkg, error_on_fail = FALSE)),
      info = paste0("does not exist (anymore): function `", pkg, "::", fn, "()`")
    )
  } else {
    warning("Package '", pkg, "' does not exist anymore")
  }
}

if (AMR:::pkg_is_available("cli")) {
  expect_true(!is.null(cli::symbol$info))
}
if (AMR:::pkg_is_available("cli")) {
}
if (AMR:::pkg_is_available("cli")) {
  expect_true(!is.null(cli::symbol$ellipsis))
}

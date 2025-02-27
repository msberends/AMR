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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# Check if these functions still exist in their package (all are in Suggests field)
# Since GitHub Actions runs every night, we will get emailed when a dependency fails based on this unit test

# functions used by import_fn()
import_functions <- c(
  "%chin%" = "data.table",
  "anti_join" = "dplyr",
  "as.data.table" = "data.table",
  "as_tibble" = "tibble",
  "chmatch" = "data.table",
  "cli_abort" = "cli",
  "cur_column" = "dplyr",
  "cur_group" = "dplyr",
  "document_position" = "rstudioapi",
  "document_range" = "rstudioapi",
  "full_join" = "dplyr",
  "getActiveDocumentContext" = "rstudioapi",
  "has_color" = "crayon",
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
  "modifyRange" = "rstudioapi",
  "new_pillar_shaft_simple" = "pillar",
  "progress_bar" = "progress",
  "read_html" = "xml2",
  "right_join" = "dplyr",
  "semi_join" = "dplyr",
  "showQuestion" = "rstudioapi",
  "symbol" = "cli",
  "tibble" = "tibble",
  "write.xlsx" = "openxlsx"
)

# functions that are called directly with ::
call_functions <- c(
  # cleaner
  "freq" = "cleaner",
  "freq.default" = "cleaner",
  "percentage" = "cleaner",
  # cli
  "symbol" = "cli",
  # curl
  "has_internet" = "curl",
  # ggplot2
  "aes" = "ggplot2",
  "arrow" = "ggplot2",
  "autoplot" = "ggplot2",
  "element_blank" = "ggplot2",
  "element_line" = "ggplot2",
  "element_text" = "ggplot2",
  "expand_limits" = "ggplot2",
  "facet_wrap" = "ggplot2",
  "fortify" = "ggplot2",
  "geom_col" = "ggplot2",
  "geom_errorbar" = "ggplot2",
  "geom_path" = "ggplot2",
  "geom_point" = "ggplot2",
  "geom_ribbon" = "ggplot2",
  "geom_segment" = "ggplot2",
  "geom_text" = "ggplot2",
  "ggplot" = "ggplot2",
  "labs" = "ggplot2",
  "position_dodge2" = "ggplot2",
  "position_fill" = "ggplot2",
  "scale_colour_discrete" = "ggplot2",
  "scale_discrete_manual" = "ggplot2",
  "scale_fill_discrete" = "ggplot2",
  "scale_fill_manual" = "ggplot2",
  "scale_x_discrete" = "ggplot2",
  "scale_y_continuous" = "ggplot2",
  "scale_y_discrete" = "ggplot2",
  "theme" = "ggplot2",
  "theme_minimal" = "ggplot2",
  "unit" = "ggplot2",
  "xlab" = "ggplot2",
  "ylab" = "ggplot2",
  # knitr
  "asis_output" = "knitr",
  "kable" = "knitr",
  "knit_print" = "knitr",
  "opts_chunk" = "knitr",
  # pillar
  "pillar_shaft" = "pillar",
  "tbl_format_footer" = "pillar",
  "tbl_sum" = "pillar",
  "type_sum" = "pillar",
  # readxl
  "read_excel" = "readxl",
  # rmarkdown
  "html_vignette" = "rmarkdown",
  # skimr
  "get_skimmers" = "skimr",
  "inline_hist" = "skimr",
  "sfl" = "skimr",
  # tibble
  "tibble" = "tibble",
  # vctrs
  "vec_arith" = "vctrs",
  "vec_cast" = "vctrs",
  "vec_math" = "vctrs",
  "vec_ptype2" = "vctrs",
  "vec_ptype_abbr" = "vctrs",
  "vec_ptype_full" = "vctrs"
)

import_functions <- c(import_functions, call_functions)

suggests <- strsplit(utils::packageDescription("AMR")$Suggests, "[,\n ]+")[[1]]
for (i in seq_len(length(import_functions))) {
  fn <- names(import_functions)[i]
  pkg <- unname(import_functions[i])
  expect_true(pkg %in% suggests,
    info = paste0("package `", pkg, "` is not in Suggests")
  )
  # function should exist in foreign pkg namespace
  if (AMR:::pkg_is_available(pkg,
    also_load = FALSE,
    min_version = if (pkg == "dplyr") "1.0.0" else NULL
  )) {
    expect_true(!is.null(AMR:::import_fn(name = fn, pkg = pkg, error_on_fail = FALSE)),
      info = paste0("Function does not exist (anymore): function `", pkg, "::", fn, "()`")
    )
  } else if (pkg != "rstudioapi") {
    warning("Package '", pkg, "' not available")
  }
}

if (AMR:::pkg_is_available("cli")) {
  expect_true(!is.null(cli::symbol$bullet) && is.character(cli::symbol$bullet) && length(cli::symbol$bullet) == 1)
  expect_true(!is.null(cli::symbol$ellipsis) && is.character(cli::symbol$ellipsis) && length(cli::symbol$ellipsis) == 1)
  expect_true(!is.null(cli::symbol$info) && is.character(cli::symbol$info) && length(cli::symbol$info) == 1)
  expect_true(!is.null(cli::symbol$sup_1) && is.character(cli::symbol$sup_1) && length(cli::symbol$sup_1) == 1)
}
if (AMR:::pkg_is_available("ggplot2")) {
  # the scale_*_mic() functions rely on these
  expect_true(is.function(ggplot2::scale_x_discrete()$transform))
  expect_true(is.function(ggplot2::scale_y_discrete()$transform))
  expect_true(is.function(ggplot2::scale_colour_discrete()$transform))
  expect_true(is.function(ggplot2::scale_fill_discrete()$transform))
}

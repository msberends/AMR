# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

###############
# NOTE TO SELF: could also have done this with the 'lifecycle' package, but why add another dependency for such an easy job?
###############

#' Lifecycles of functions in the `AMR` package
#' @name lifecycle
#' @rdname lifecycle
#' @description Our functions are categorised using [the lifecycle circle of the `tidyverse` as found on www.tidyverse.org/lifecycle](https://www.tidyverse.org/lifecycle).
#' 
#' \if{html}{\figure{lifecycle_tidyverse.svg}{options: height=200px style=margin-bottom:5px} \cr}
#' This page contains a section for every lifecycle (with text borrowed from the aforementioned `tidyverse` website), so they can be used in the manual pages of our functions. 
#' @section Experimental lifecycle:
#' \if{html}{\figure{lifecycle_experimental.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **experimental**. An experimental function is in the very early stages of development. The unlying code might be changing frequently as we rapidly iterate and explore variations in search of the best fit. Experimental functions might be removed without deprecation, so you are generally best off waiting until a function is more mature before you use it in production code. Experimental functions will not be included in releases we submit to CRAN.
#' @section Maturing lifecycle:
#' \if{html}{\figure{lifecycle_maturing.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **maturing**. The unlying code of a maturing function has been roughed out, but finer details might still change. We will strive to maintain backward compatibility, but the function needs wider usage and more extensive testing in order to optimise the unlying code.
#' @section Stable lifecycle:
#' \if{html}{\figure{lifecycle_stable.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **stable**. In a stable function, we are largely happy with the unlying code, and major changes are unlikely. This means that the unlying code will generally evolve by adding new arguments; we will avoid removing arguments or changing the meaning of existing arguments.
#' 
#' If the unlying code needs breaking changes, they will occur gradually. To begin with, the function or argument will be deprecated; it will continue to work but will emit an message informing you of the change. Next, typically after at least one newly released version on CRAN, the message will be transformed to an error.
#' @section Retired lifecycle:
#' \if{html}{\figure{lifecycle_retired.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **retired**. A retired function is no longer under active development, and (if appropiate) a better alternative is available. We will only make the necessary changes to ensure that retired functions remain available. No new arguments will be added, and only the most critical bugs will be fixed.
#' @section Archived lifecycle:
#' \if{html}{\figure{lifecycle_archived.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **archived**. The development of an archived function has ended, and it is no longer available in future package versions.
#' @section Dormant lifecycle:
#' \if{html}{\figure{lifecycle_dormant.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **dormant**. A dormant function is currently not under active development and has not reached a stable phase. We might return to it in the future. As with experimental functions, you are best off waiting until a function is more mature before you use it in production code.
#' @section Questioning lifecycle:
#' \if{html}{\figure{lifecycle_questioning.svg}{options: style=margin-bottom:5px} \cr}
#' The [lifecycle][AMR::lifecycle] of this function is **questioning**. We are no longer convinced that this function is the optimal approach (but we do not know yet what a better approach would be), or whether this function should be in our `AMR` package at all.
NULL

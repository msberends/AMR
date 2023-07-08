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

expect_true(as.mic(8) == as.mic("8"))
expect_true(as.mic("1") > as.mic("<=0.0625"))
expect_true(as.mic("1") < as.mic(">=32"))
expect_true(is.mic(as.mic(8)))
# expect_true(as.mic(1024) < as.mic(">1024"))
# expect_true(as.mic("<1024") > as.mic("1024"))


expect_equal(as.double(as.mic(">=32")), 32)
expect_equal(as.numeric(as.mic(">=32")), 32)
expect_equal(
  as.integer(as.mic(">=32")), # should be factor level, not the MIC
  as.integer(factor(as.character(">=32"),
    levels = levels(as.mic(">=32"))
  ))
)
expect_equal(suppressWarnings(as.logical(as.mic("INVALID VALUE"))), NA)

# all levels should be valid MICs
x <- as.mic(c(2, 4))
expect_inherits(x[1], "mic")
expect_inherits(x[[1]], "mic")
expect_inherits(c(x[1], x[9]), "mic")
expect_inherits(unique(x[1], x[9]), "mic")
expect_inherits(droplevels(c(x[1], x[9]), as.mic = TRUE), "factor")
expect_inherits(droplevels(c(x[1], x[9]), as.mic = TRUE), "mic")
x[2] <- 32
expect_inherits(x, "mic")
# expect_warning(as.mic("INVALID VALUE"))

pdf(NULL) # prevent Rplots.pdf being created
expect_silent(barplot(as.mic(c(1, 2, 4, 8))))
expect_silent(plot(as.mic(c(1, 2, 4, 8))))
expect_silent(plot(as.mic(c(1, 2, 4, 8)), expand = FALSE))
expect_silent(plot(as.mic(c(1, 2, 4, 8)), mo = "Escherichia coli", ab = "cipr"))
if (AMR:::pkg_is_available("ggplot2")) {
  expect_inherits(ggplot2::autoplot(as.mic(c(1, 2, 4, 8))), "gg")
  expect_inherits(ggplot2::autoplot(as.mic(c(1, 2, 4, 8)), expand = FALSE), "gg")
  expect_inherits(ggplot2::autoplot(as.mic(c(1, 2, 4, 8, 32)), mo = "Escherichia coli", ab = "cipr"), "gg")
}
expect_stdout(print(as.mic(c(1, 2, 4, 8))))

expect_inherits(summary(as.mic(c(2, 8))), c("summaryDefault", "table"))

if (AMR:::pkg_is_available("tibble")) {
  expect_stdout(print(tibble::tibble(m = as.mic(2:4))))
}

# all mathematical operations
x <- random_mic(50)
x_double <- as.double(gsub("[<=>]+", "", as.character(x)))
suppressWarnings(expect_identical(mean(x), mean(x_double)))
suppressWarnings(expect_identical(median(x), median(x_double)))
suppressWarnings(expect_identical(quantile(x), quantile(x_double)))
suppressWarnings(expect_identical(abs(x), abs(x_double)))
suppressWarnings(expect_identical(sign(x), sign(x_double)))
suppressWarnings(expect_identical(sqrt(x), sqrt(x_double)))
suppressWarnings(expect_identical(floor(x), floor(x_double)))
suppressWarnings(expect_identical(ceiling(x), ceiling(x_double)))
suppressWarnings(expect_identical(trunc(x), trunc(x_double)))
suppressWarnings(expect_identical(round(x), round(x_double)))
suppressWarnings(expect_identical(signif(x), signif(x_double)))
suppressWarnings(expect_identical(exp(x), exp(x_double)))
suppressWarnings(expect_identical(log(x), log(x_double)))
suppressWarnings(expect_identical(log10(x), log10(x_double)))
suppressWarnings(expect_identical(log2(x), log2(x_double)))
suppressWarnings(expect_identical(expm1(x), expm1(x_double)))
suppressWarnings(expect_identical(log1p(x), log1p(x_double)))
suppressWarnings(expect_identical(cos(x), cos(x_double)))
suppressWarnings(expect_identical(sin(x), sin(x_double)))
suppressWarnings(expect_identical(tan(x), tan(x_double)))
if (getRversion() >= "3.1") {
  suppressWarnings(expect_identical(cospi(x), cospi(x_double)))
  suppressWarnings(expect_identical(sinpi(x), sinpi(x_double)))
  suppressWarnings(expect_identical(tanpi(x), tanpi(x_double)))
}
suppressWarnings(expect_identical(acos(x), acos(x_double)))
suppressWarnings(expect_identical(asin(x), asin(x_double)))
suppressWarnings(expect_identical(atan(x), atan(x_double)))
suppressWarnings(expect_identical(cosh(x), cosh(x_double)))
suppressWarnings(expect_identical(sinh(x), sinh(x_double)))
suppressWarnings(expect_identical(tanh(x), tanh(x_double)))
suppressWarnings(expect_identical(acosh(x), acosh(x_double)))
suppressWarnings(expect_identical(asinh(x), asinh(x_double)))
suppressWarnings(expect_identical(atanh(x), atanh(x_double)))
suppressWarnings(expect_identical(lgamma(x), lgamma(x_double)))
suppressWarnings(expect_identical(gamma(x), gamma(x_double)))
suppressWarnings(expect_identical(digamma(x), digamma(x_double)))
suppressWarnings(expect_identical(trigamma(x), trigamma(x_double)))
suppressWarnings(expect_identical(cumsum(x), cumsum(x_double)))
suppressWarnings(expect_identical(cumprod(x), cumprod(x_double)))
suppressWarnings(expect_identical(cummax(x), cummax(x_double)))
suppressWarnings(expect_identical(cummin(x), cummin(x_double)))
suppressWarnings(expect_identical(!x, !x_double))

suppressWarnings(expect_identical(all(x), all(x_double)))
suppressWarnings(expect_identical(any(x), any(x_double)))
suppressWarnings(expect_identical(sum(x), sum(x_double)))
suppressWarnings(expect_identical(prod(x), prod(x_double)))
suppressWarnings(expect_identical(min(x), min(x_double)))
suppressWarnings(expect_identical(max(x), max(x_double)))
suppressWarnings(expect_identical(range(x), range(x_double)))

el1 <- random_mic(50)
el1_double <- as.double(gsub("[<=>]+", "", as.character(el1)))
el2 <- random_mic(50)
el2_double <- as.double(gsub("[<=>]+", "", as.character(el2)))
suppressWarnings(expect_identical(el1 + el2, el1_double + el2_double))
suppressWarnings(expect_identical(el1 - el2, el1_double - el2_double))
suppressWarnings(expect_identical(el1 * el2, el1_double * el2_double))
suppressWarnings(expect_identical(el1 / el2, el1_double / el2_double))
suppressWarnings(expect_identical(el1^el2, el1_double^el2_double))
suppressWarnings(expect_identical(el1 %% el2, el1_double %% el2_double))
suppressWarnings(expect_identical(el1 %/% el2, el1_double %/% el2_double))
suppressWarnings(expect_identical(el1 & el2, el1_double & el2_double))
suppressWarnings(expect_identical(el1 | el2, el1_double | el2_double))
suppressWarnings(expect_identical(el1 == el2, el1_double == el2_double))
suppressWarnings(expect_identical(el1 != el2, el1_double != el2_double))
suppressWarnings(expect_identical(el1 < el2, el1_double < el2_double))
suppressWarnings(expect_identical(el1 <= el2, el1_double <= el2_double))
suppressWarnings(expect_identical(el1 >= el2, el1_double >= el2_double))
suppressWarnings(expect_identical(el1 > el2, el1_double > el2_double))

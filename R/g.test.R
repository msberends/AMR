# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' \emph{G}-test of matrix or vector
#'
#' A \emph{G}-test can be used to see whether the number of observations in each category fits a theoretical expectation (called a \strong{\emph{G}-test of goodness-of-fit}), or to see whether the proportions of one variable are different for different values of the other variable (called a \strong{\emph{G}-test of independence}).
#' @param x numeric vector or matrix
#' @param y expected value of \code{x}. Leave empty to determine automatically. This can also be ratios of \code{x}, e.g. calculated with \code{\link{vector2ratio}}.
#' @param alpha value to test the p value against
#' @param info logical to determine whether the analysis should be printed
#' @param minimum the test with fail if any of the observed values is below this value. Use \code{minimum = 30} for microbial epidemiology, to prevent calculating a p value when less than 30 isolates are available.
#' @section \emph{G}-test of goodness-of-fit (likelihood ratio test):
#' Use the \emph{G}-test of goodness-of-fit when you have one nominal variable with two or more values (such as male and female, or red, pink and white flowers). You compare the observed counts of numbers of observations in each category with the expected counts, which you calculate using some kind of theoretical expectation (such as a 1:1 sex ratio or a 1:2:1 ratio in a genetic cross).
#'
#' If the expected number of observations in any category is too small, the \emph{G}-test may give inaccurate results, and you should use an exact test instead. See the web page on small sample sizes for discussion of what "small" means.
#'
#' The \emph{G}-test of goodness-of-fit is an alternative to the chi-square test of goodness-of-fit; each of these tests has some advantages and some disadvantages, and the results of the two tests are usually very similar.
#'
#' @section \emph{G}-test of independence:
#' Use the \emph{G}-test of independence when you have two nominal variables, each with two or more possible values. You want to know whether the proportions for one variable are different among values of the other variable.
#'
#' It is also possible to do a \emph{G}-test of independence with more than two nominal variables. For example, Jackson et al. (2013) also had data for children under 3, so you could do an analysis of old vs. young, thigh vs. arm, and reaction vs. no reaction, all analyzed together.
#'
#' Fisher's exact test is more accurate than the \emph{G}-test of independence when the expected numbers are small, so it is recommend to only use the \emph{G}-test if your total sample size is greater than 1000.
#'
#' The \emph{G}-test of independence is an alternative to the chi-square test of independence, and they will give approximately the same results.
#' @section How the test works:
#' Unlike the exact test of goodness-of-fit, the \emph{G}-test does not directly calculate the probability of obtaining the observed results or something more extreme. Instead, like almost all statistical tests, the \emph{G}-test has an intermediate step; it uses the data to calculate a test statistic that measures how far the observed data are from the null expectation. You then use a mathematical relationship, in this case the chi-square distribution, to estimate the probability of obtaining that value of the test statistic.
#'
#' The \emph{G}-test uses the log of the ratio of two likelihoods as the test statistic, which is why it is also called a likelihood ratio test or log-likelihood ratio test. The formula to calculate a \emph{G}-statistic is:
#'
#' \code{G <- 2 * sum(x * log(x / x.expected))}
#'
#' Since this is chi-square distributed, the p value can be calculated with:
#'
#' \code{p <- 1 - stats::pchisq(G, df))}
#'
#' where \code{df} are the degrees of freedom: \code{max(NROW(x) - 1, 1) * max(NCOL(x) - 1, 1)}.
#'
#' If there are more than two categories and you want to find out which ones are significantly different from their null expectation, you can use the same method of testing each category vs. the sum of all categories, with the Bonferroni correction. You use \emph{G}-tests for each category, of course.
#' @keywords chi
#' @seealso \code{\link{chisq.test}}
#' @references McDonald, J.H. 2014. \strong{Handbook of Biological Statistics (3rd ed.)}. Sparky House Publishing, Baltimore, Maryland. \url{http://www.biostathandbook.com/gtestgof.html}.
#' @export
#' @importFrom stats pchisq
#' @importFrom dplyr %>%
#' @examples
#' # = EXAMPLE 1 =
#' # Shivrain et al. (2006) crossed clearfield rice (which are resistant
#' # to the herbicide imazethapyr) with red rice (which are susceptible to
#' # imazethapyr). They then crossed the hybrid offspring and examined the
#' # F2 generation, where they found 772 resistant plants, 1611 moderately
#' # resistant plants, and 737 susceptible plants. If resistance is controlled
#' # by a single gene with two co-dominant alleles, you would expect a 1:2:1
#' # ratio.
#'
#' x <- c(772, 1611, 737)
#' x.expected <- vector2ratio(x, ratio = "1:2:1")
#' x.expected
#' # 780 1560 780
#'
#' g.test(x, x.expected)
#' # p = 0.12574.
#'
#' # There is no significant difference from a 1:2:1 ratio.
#' # Meaning: resistance controlled by a single gene with two co-dominant
#' # alleles, is plausible.
#'
#'
#' # = EXAMPLE 2 =
#' # Red crossbills (Loxia curvirostra) have the tip of the upper bill either
#' # right or left of the lower bill, which helps them extract seeds from pine
#' # cones. Some have hypothesized that frequency-dependent selection would
#' # keep the number of right and left-billed birds at a 1:1 ratio. Groth (1992)
#' # observed 1752 right-billed and 1895 left-billed crossbills.
#'
#' x <- c(1752, 1895)
#' x.expected <- vector2ratio(x, ratio = c(1, 1))
#' x.expected
#' # 1823.5 1823.5
#'
#' g.test(x, x.expected)
#' # p = 0.01787343
#'
#' # There is a significant difference from a 1:1 ratio.
#' # Meaning: there are significantly more left-billed birds.
#'
g.test <- function(x,
                   y = NULL,
                   alpha = 0.05,
                   info = TRUE,
                   minimum = 0) {

  if (sum(x) < 1000) {
    warning('the sum of all observations is < 1000, consider using a Fishers Exact test instead.')
  }

  if (!is.numeric(x)) {
    stop('`x` must be a vector or matrix with numeric values.')
  }
  if (!is.matrix(x) & is.null(y)) {
    stop('if `x` is not a matrix, `y` must be given.')
  }

  # if (!is.matrix(x)) {
  #   x <- matrix(x, dimnames = list(rep("", length(x)), ""))
  # }

  x.expected <- y
  # if (!is.null(x.expected) & !is.matrix(x.expected)) {
  #   x.expected <- matrix(x.expected, dimnames = list(rep("", length(x.expected)), ""))
  # }

  if (NCOL(x) > 1) {
    matrix2tbl <- function(x) {
      if (!is.matrix(x)) {
        x <- matrix(x, dimnames = list(rep("", length(x)), ""))
      }
      x <- rbind(cbind(x, rowSums(x)), colSums(cbind(x, rowSums(x)))) %>%
        as.data.frame()
      colnames(x) <- c(paste0('c', 1:(NCOL(x) - 1)), 'Total')
      rownames(x) <- c(strrep(" ", 1:(NROW(x) - 1)), 'Total')
      x
    }

    x.with_totals <- matrix2tbl(x)
    if (is.null(x.expected)) {
      x.expected <- outer(rowSums(x), colSums(x), "*") / sum(x)
      x.expected.with_totals <- matrix2tbl(x.expected)
    }

  } else {
    x.with_totals <- x
    if (is.null(x.expected)) {
      x.expected <- matrix(rep(mean(x), length(x)), dimnames = list(rep("", length(x)), ""))
      warning('Expected values set to the mean of `x`, because it consists of only one column.')
    }
    x.expected.with_totals <- x.expected
  }

  if (any(x < minimum)) {
    warning('One of the observed values is lower than the required minimum of ', minimum, '.')
    return(NA_real_)
  }
  if (any(x.expected < 5)) {
    warning('One of the expected values is lower than 5, thus G-statistic approximation may be incorrect.',
            '\n  Consider doing an Exact test (exact.test).')
  }

  Gstat <- 2 * base::sum(x * log(x / x.expected), na.rm = TRUE)
  df <- base::max(NROW(x) - 1, 1, na.rm = TRUE) * base::max(NCOL(x) - 1, 1, na.rm = TRUE)

  pval <- 1 - stats::pchisq(Gstat, df = df)

  if (info == TRUE) {

    if (NROW(x) == 2 & NCOL(x) == 2) {
      cat('G-test of independence\n\n')
    } else {
      cat('G-test of goodness-of-fit\n')
      cat('(likelihood ratio test)\n\n')
    }

    cat(paste0("(O) Observed values:\n"))
    print(x.with_totals %>% round(2))

    cat('\n')
    cat('(E) Expected values under null hypothesis:\n')
    print(x.expected.with_totals %>% round(2))

    cat('\n============[G-test]============\n')
    cat(' G-statistic         :', round(Gstat, 4), '\n')
    cat(' Degrees of freedom  :', df, '\n')
    cat(' P value             :', round(pval, 4), '\n')
    cat(' Alpha               :', round(alpha, 2), '\n')
    cat(' Significance        :', p.symbol(pval, "-"), '\n')
    cat('================================\n\n')

  }

  pval

}

#' Transform vector to ratio
#' @param x vector of values
#' @param ratio vector with ratios of \code{x} and with same length (like \code{ratio = c(1, 2, 1)}) or a text with characters \code{":"}, \code{"-"} or \code{","} (like \code{ratio = "1:2:1"} or even \code{ratio = "1:2:1.25"})
#' @export
#' @seealso \code{\link{g.test}}
#' @references McDonald, J.H. 2014. \strong{Handbook of Biological Statistics (3rd ed.)}. Sparky House Publishing, Baltimore, Maryland.
#' @importFrom dplyr %>%
#' @examples
#' # = EXAMPLE 1 =
#' # Shivrain et al. (2006) crossed clearfield rice (which are resistant
#' # to the herbicide imazethapyr) with red rice (which are susceptible to
#' # imazethapyr). They then crossed the hybrid offspring and examined the
#' # F2 generation, where they found 772 resistant plants, 1611 moderately
#' # resistant plants, and 737 susceptible plants. If resistance is controlled
#' # by a single gene with two co-dominant alleles, you would expect a 1:2:1
#' # ratio.
#'
#' x <- c(772, 1611, 737)
#' x.expected <- vector2ratio(x, ratio = "1:2:1")
#' x.expected
#' # 780 1560 780
#'
#' g.test(x, x.expected)
#' # p = 0.12574.
#'
#' # There is no significant difference from a 1:2:1 ratio.
#' # Meaning: resistance controlled by a single gene with two co-dominant
#' # alleles, is plausible.
#'
#'
#' # = EXAMPLE 2 =
#' # Red crossbills (Loxia curvirostra) have the tip of the upper bill either
#' # right or left of the lower bill, which helps them extract seeds from pine
#' # cones. Some have hypothesized that frequency-dependent selection would
#' # keep the number of right and left-billed birds at a 1:1 ratio. Groth (1992)
#' # observed 1752 right-billed and 1895 left-billed crossbills.
#'
#' x <- c(1752, 1895)
#' x.expected <- vector2ratio(x, ratio = c(1, 1))
#' x.expected
#' # 1823.5 1823.5
#'
#' g.test(x, x.expected)
#' # p = 0.01787343
#'
#' # There is a significant difference from a 1:1 ratio.
#' # Meaning: there are significantly more left-billed birds.
#'
vector2ratio <- function(x, ratio) {
  if (!all(is.numeric(x))) {
    stop('`x` must be a vector of numeric values.')
  }
  if (length(ratio) == 1) {
    if (ratio %like% '^([0-9]+([.][0-9]+)?[-,:])+[0-9]+([.][0-9]+)?$') {
      # support for "1:2:1", "1-2-1", "1,2,1" and even "1.75:2:1.5"
      ratio <- ratio %>% base::strsplit("[-,:]") %>% base::unlist() %>% base::as.double()
    } else {
      stop('Invalid `ratio`: ', ratio, '.')
    }
  }
  if (length(x) != length(ratio)) {
    stop('`x` and `ratio` must be of same size.')
  }
  base::sum(x, na.rm = TRUE) * (ratio / base::sum(ratio, na.rm = TRUE))
}

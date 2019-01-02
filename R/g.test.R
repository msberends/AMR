# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' \emph{G}-test for Count Data
#'
#' \code{g.test} performs chi-squared contingency table tests and goodness-of-fit tests, just like \code{\link{chisq.test}} but is more reliable [1]. A \emph{G}-test can be used to see whether the number of observations in each category fits a theoretical expectation (called a \strong{\emph{G}-test of goodness-of-fit}), or to see whether the proportions of one variable are different for different values of the other variable (called a \strong{\emph{G}-test of independence}).
#' @inherit stats::chisq.test params return
#' @details If \code{x} is a matrix with one row or column, or if \code{x} is a vector and \code{y} is not given, then a \emph{goodness-of-fit test} is performed (\code{x} is treated as a one-dimensional contingency table). The entries of \code{x} must be non-negative integers. In this case, the hypothesis tested is whether the population probabilities equal those in \code{p}, or are all equal if \code{p} is not given.
#'
#'   If \code{x} is a matrix with at least two rows and columns, it is taken as a two-dimensional contingency table: the entries of \code{x} must be non-negative integers.  Otherwise, \code{x} and \code{y} must be vectors or factors of the same length; cases with missing values are removed, the objects are coerced to factors, and the contingency table is computed from these.  Then Pearson's chi-squared test is performed of the null hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency table is the product of the row and column marginals.
#'
#'   The p-value is computed from the asymptotic chi-squared distribution of the test statistic.
#'
#'   In the contingency table case simulation is done by random sampling from the set of all contingency tables with given marginals, and works only if the marginals are strictly positive. Note that this is not the usual sampling situation assumed for a chi-squared test (like the \emph{G}-test) but rather that for Fisher's exact test.
#'
#'   In the goodness-of-fit case simulation is done by random sampling from the discrete distribution specified by \code{p}, each sample being of size \code{n = sum(x)}. This simulation is done in \R and may be slow.
#' @section \emph{G}-test of goodness-of-fit (likelihood ratio test):
#' Use the \emph{G}-test of goodness-of-fit when you have one nominal variable with two or more values (such as male and female, or red, pink and white flowers). You compare the observed counts of numbers of observations in each category with the expected counts, which you calculate using some kind of theoretical expectation (such as a 1:1 sex ratio or a 1:2:1 ratio in a genetic cross).
#'
#' If the expected number of observations in any category is too small, the \emph{G}-test may give inaccurate results, and you should use an exact test instead (\code{\link{fisher.test}}).
#'
#' The \emph{G}-test of goodness-of-fit is an alternative to the chi-square test of goodness-of-fit (\code{\link{chisq.test}}); each of these tests has some advantages and some disadvantages, and the results of the two tests are usually very similar.
#'
#' @section \emph{G}-test of independence:
#' Use the \emph{G}-test of independence when you have two nominal variables, each with two or more possible values. You want to know whether the proportions for one variable are different among values of the other variable.
#'
#' It is also possible to do a \emph{G}-test of independence with more than two nominal variables. For example, Jackson et al. (2013) also had data for children under 3, so you could do an analysis of old vs. young, thigh vs. arm, and reaction vs. no reaction, all analyzed together.
#'
#' Fisher's exact test (\code{\link{fisher.test}}) is more accurate than the \emph{G}-test of independence when the expected numbers are small, so it is recommend to only use the \emph{G}-test if your total sample size is greater than 1000.
#'
#' The \emph{G}-test of independence is an alternative to the chi-square test of independence (\code{\link{chisq.test}}), and they will give approximately the same results.
#' @section How the test works:
#' Unlike the exact test of goodness-of-fit (\code{\link{fisher.test}}), the \emph{G}-test does not directly calculate the probability of obtaining the observed results or something more extreme. Instead, like almost all statistical tests, the \emph{G}-test has an intermediate step; it uses the data to calculate a test statistic that measures how far the observed data are from the null expectation. You then use a mathematical relationship, in this case the chi-square distribution, to estimate the probability of obtaining that value of the test statistic.
#'
#' The \emph{G}-test uses the log of the ratio of two likelihoods as the test statistic, which is why it is also called a likelihood ratio test or log-likelihood ratio test. The formula to calculate a \emph{G}-statistic is:
#'
#' \code{G <- 2 * sum(x * log(x / E))}
#'
#' where \code{E} are the expected values. Since this is chi-square distributed, the p value can be calculated with:
#'
#' \code{p <- stats::pchisq(G, df, lower.tail = FALSE)}
#'
#' where \code{df} are the degrees of freedom.
#'
#' If there are more than two categories and you want to find out which ones are significantly different from their null expectation, you can use the same method of testing each category vs. the sum of all categories, with the Bonferroni correction. You use \emph{G}-tests for each category, of course.
#' @keywords chi
#' @seealso \code{\link{chisq.test}}
#' @references [1] McDonald, J.H. 2014. \strong{Handbook of Biological Statistics (3rd ed.)}. Sparky House Publishing, Baltimore, Maryland. \url{http://www.biostathandbook.com/gtestgof.html}.
#' @source This code is almost identical to \code{\link{chisq.test}}, except that:
#' \itemize{
#'   \item{The calculation of the statistic was changed to \code{2 * sum(x * log(x / E))}}
#'   \item{Yates' continuity correction was removed as it does not apply to a \emph{G}-test}
#'   \item{The possibility to simulate p values with \code{simulate.p.value} was removed}
#' }
#' @export
#' @importFrom stats pchisq complete.cases
#' @inheritSection AMR Read more on our website!
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
#' G <- g.test(x, p = c(1, 2, 1) / 4)
#' # G$p.value = 0.12574.
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
#' g.test(x)
#' # p = 0.01787343
#'
#' # There is a significant difference from a 1:1 ratio.
#' # Meaning: there are significantly more left-billed birds.
#'
g.test <- function(x,
                    y = NULL,
                    # correct = TRUE,
                    p = rep(1/length(x), length(x)),
                    rescale.p = FALSE) {
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      x <- as.vector(x)
  }
  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME2 <- deparse(substitute(y))
    xname <- if (length(DNAME) > 1L || nchar(DNAME, "w") >
                 30)
      ""
    else DNAME
    yname <- if (length(DNAME2) > 1L || nchar(DNAME2, "w") >
                 30)
      ""
    else DNAME2
    OK <- complete.cases(x, y)
    x <- factor(x[OK])
    y <- factor(y[OK])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L))
      stop("'x' and 'y' must have at least 2 levels")
    x <- table(x, y)
    names(dimnames(x)) <- c(xname, yname)
    DNAME <- paste(paste(DNAME, collapse = "\n"), "and",
                   paste(DNAME2, collapse = "\n"))
  }
  if (any(x < 0) || anyNA(x))
    stop("all entries of 'x' must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of 'x' must be positive")
  # if (simulate.p.value) {
  #   setMETH <- function() METHOD <<- paste(METHOD, "with simulated p-value\n\t (based on",
  #                                          B, "replicates)")
  #   almost.1 <- 1 - 64 * .Machine$double.eps
  # }
  if (is.matrix(x)) {
    METHOD <- "G-test of independence"
    nr <- as.integer(nrow(x))
    nc <- as.integer(ncol(x))
    if (is.na(nr) || is.na(nc) || is.na(nr * nc))
      stop("invalid nrow(x) or ncol(x)", domain = NA)
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc, "*")/n
    v <- function(r, c, n) c * r * (n - r) * (n - c)/n^3
    V <- outer(sr, sc, v, n)
    dimnames(E) <- dimnames(x)
    # if (simulate.p.value && all(sr > 0) && all(sc > 0)) {
    #   setMETH()
    #   tmp <- .Call(chisq_sim, sr, sc, B, E, PACKAGE = "stats")
    #   STATISTIC <- 2 * sum(x * log(x / E)) # sum(sort((x - E)^2/E, decreasing = TRUE)) for chisq.test
    #   PARAMETER <- NA
    #   PVAL <- (1 + sum(tmp >= almost.1 * STATISTIC))/(B +
    #                                                     1)
    # }
    # else {
      # if (simulate.p.value)
      #   warning("cannot compute simulated p-value with zero marginals")
      # if (correct && nrow(x) == 2L && ncol(x) == 2L) {
      #   YATES <- min(0.5, abs(x - E))
      #   if (YATES > 0)
      #     METHOD <- paste(METHOD, "with Yates' continuity correction")
      # }
      # else YATES <- 0
      STATISTIC <- 2 * sum(x * log(x / E)) # sum((abs(x - E) - YATES)^2/E) for chisq.test
      PARAMETER <- (nr - 1L) * (nc - 1L)
      PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
   # }
  }
  else {
    if (length(dim(x)) > 2L)
      stop("invalid 'x'")
    if (length(x) == 1L)
      stop("'x' must at least have 2 elements")
    if (length(x) != length(p))
      stop("'x' and 'p' must have the same number of elements")
    if (any(p < 0))
      stop("probabilities must be non-negative.")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
      if (rescale.p)
        p <- p/sum(p)
      else stop("probabilities must sum to 1.")
    }
    METHOD <- "G-test of goodness-of-fit (likelihood ratio test)"
    E <- n * p
    V <- n * p * (1 - p)
    STATISTIC <- 2 * sum(x * log(x / E)) # sum((x - E)^2/E) for chisq.test
    names(E) <- names(x)
    # if (simulate.p.value) {
    #   setMETH()
    #   nx <- length(x)
    #   sm <- matrix(sample.int(nx, B * n, TRUE, prob = p),
    #                nrow = n)
    #   ss <- apply(sm, 2L, function(x, E, k) {
    #     sum((table(factor(x, levels = 1L:k)) - E)^2/E)
    #   }, E = E, k = nx)
    #   PARAMETER <- NA
    #   PVAL <- (1 + sum(ss >= almost.1 * STATISTIC))/(B +
    #                                                    1)
    # }
    # else {
      PARAMETER <- length(x) - 1
      PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    # }
  }
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  # if (any(E < 5) && is.finite(PARAMETER))
  #   warning("G-statistic approximation may be incorrect")

  # suggest fisher.test when total is < 1000 (John McDonald, Handbook of Biological Statistics, 2014)
  if (sum(x, na.rm = TRUE) < 1000 && is.finite(PARAMETER)) {
    warning("G-statistic approximation may be incorrect, consider Fisher's Exact test")
  } else if (any(E < 5) && is.finite(PARAMETER)) {
    warning("G-statistic approximation may be incorrect, consider Fisher's Exact test")
  }
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME,
                 observed = x, expected = E, residuals = (x - E)/sqrt(E),
                 stdres = (x - E)/sqrt(V)), class = "htest")
}

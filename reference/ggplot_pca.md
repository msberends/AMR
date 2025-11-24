# PCA Biplot with `ggplot2`

Produces a `ggplot2` variant of a so-called
[biplot](https://en.wikipedia.org/wiki/Biplot) for PCA (principal
component analysis), but is more flexible and more appealing than the
base R [`biplot()`](https://rdrr.io/r/stats/biplot.html) function.

## Usage

``` r
ggplot_pca(x, choices = 1:2, scale = 1, pc.biplot = TRUE,
  labels = NULL, labels_textsize = 3, labels_text_placement = 1.5,
  groups = NULL, ellipse = TRUE, ellipse_prob = 0.68,
  ellipse_size = 0.5, ellipse_alpha = 0.5, points_size = 2,
  points_alpha = 0.25, arrows = TRUE, arrows_colour = "darkblue",
  arrows_size = 0.5, arrows_textsize = 3, arrows_textangled = TRUE,
  arrows_alpha = 0.75, base_textsize = 10, ...)
```

## Source

The `ggplot_pca()` function is based on the `ggbiplot()` function from
the `ggbiplot` package by Vince Vu, as found on GitHub:
<https://github.com/vqv/ggbiplot> (retrieved: 2 March 2020, their latest
commit:
[`7325e88`](https://github.com/vqv/ggbiplot/commit/7325e880485bea4c07465a0304c470608fffb5d9);
12 February 2015).

As per their GPL-2 licence that demands documentation of code changes,
the changes made based on the source code were:

1.  Rewritten code to remove the dependency on packages `plyr`, `scales`
    and `grid`

2.  Parametrised more options, like arrow and ellipse settings

3.  Hardened all input possibilities by defining the exact type of user
    input for every argument

4.  Added total amount of explained variance as a caption in the plot

5.  Cleaned all syntax based on the `lintr` package, fixed grammatical
    errors and added integrity checks

6.  Updated documentation

## Arguments

- x:

  An object returned by
  [`pca()`](https://amr-for-r.org/reference/pca.md),
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) or
  [`princomp()`](https://rdrr.io/r/stats/princomp.html).

- choices:

  length 2 vector specifying the components to plot. Only the default is
  a biplot in the strict sense.

- scale:

  The variables are scaled by `lambda ^ scale` and the observations are
  scaled by `lambda ^ (1-scale)` where `lambda` are the singular values
  as computed by [`princomp`](https://rdrr.io/r/stats/princomp.html).
  Normally `0 <= scale <= 1`, and a warning will be issued if the
  specified `scale` is outside this range.

- pc.biplot:

  If true, use what Gabriel (1971) refers to as a "principal component
  biplot", with `lambda = 1` and observations scaled up by sqrt(n) and
  variables scaled down by sqrt(n). Then inner products between
  variables approximate covariances and distances between observations
  approximate Mahalanobis distance.

- labels:

  An optional vector of labels for the observations. If set, the labels
  will be placed below their respective points. When using the
  [`pca()`](https://amr-for-r.org/reference/pca.md) function as input
  for `x`, this will be determined automatically based on the attribute
  `non_numeric_cols`, see
  [`pca()`](https://amr-for-r.org/reference/pca.md).

- labels_textsize:

  The size of the text used for the labels.

- labels_text_placement:

  Adjustment factor the placement of the variable names (`>=1` means
  further away from the arrow head).

- groups:

  An optional vector of groups for the labels, with the same length as
  `labels`. If set, the points and labels will be coloured according to
  these groups. When using the
  [`pca()`](https://amr-for-r.org/reference/pca.md) function as input
  for `x`, this will be determined automatically based on the attribute
  `non_numeric_cols`, see
  [`pca()`](https://amr-for-r.org/reference/pca.md).

- ellipse:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  normal data ellipse should be drawn for each group (set with
  `groups`).

- ellipse_prob:

  Statistical size of the ellipse in normal probability.

- ellipse_size:

  The size of the ellipse line.

- ellipse_alpha:

  The alpha (transparency) of the ellipse line.

- points_size:

  The size of the points.

- points_alpha:

  The alpha (transparency) of the points.

- arrows:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  arrows should be drawn.

- arrows_colour:

  The colour of the arrow and their text.

- arrows_size:

  The size (thickness) of the arrow lines.

- arrows_textsize:

  The size of the text at the end of the arrows.

- arrows_textangled:

  A [logical](https://rdrr.io/r/base/logical.html) whether the text at
  the end of the arrows should be angled.

- arrows_alpha:

  The alpha (transparency) of the arrows and their text.

- base_textsize:

  The text size for all plot elements except the labels and arrows.

- ...:

  Arguments passed on to functions.

## Details

The colours for labels and points can be changed by adding another scale
layer for colour, such as
[`scale_colour_viridis_d()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
and
[`scale_colour_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html).

## Examples

``` r
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.

# \donttest{
if (require("dplyr")) {
  # calculate the resistance per group first
  resistance_data <- example_isolates %>%
    group_by(
      order = mo_order(mo), # group on anything, like order
      genus = mo_genus(mo)
    ) %>% #   and genus as we do here;
    filter(n() >= 30) %>% # filter on only 30 results per group
    summarise_if(is.sir, resistance) # then get resistance of all drugs

  # now conduct PCA for certain antimicrobial drugs
  pca_result <- resistance_data %>%
    pca(AMC, CXM, CTX, CAZ, GEN, TOB, TMP, SXT)

  summary(pca_result)

  # old base R plotting method:
  biplot(pca_result, main = "Base R biplot")

  # new ggplot2 plotting method using this package:
  if (require("ggplot2")) {
    ggplot_pca(pca_result) +
      labs(title = "ggplot2 biplot")
  }
  if (require("ggplot2")) {
    # still extendible with any ggplot2 function
    ggplot_pca(pca_result) +
      scale_colour_viridis_d() +
      labs(title = "ggplot2 biplot")
  }
}
#> Warning: There were 73 warnings in `summarise()`.
#> The first warning was:
#> ℹ In argument: `PEN = (function (..., minimum = 30, as_percent = FALSE,
#>   only_all_tested = FALSE) ...`.
#> ℹ In group 5: `order = "Lactobacillales"` `genus = "Enterococcus"`.
#> Caused by warning:
#> ! Introducing NA: only 14 results available for PEN in group: order =
#> "Lactobacillales", genus = "Enterococcus" (`minimum` = 30).
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 72 remaining warnings.
#> ℹ Columns selected for PCA: "AMC", "CAZ", "CTX", "CXM", "GEN", "SXT",
#>   "TMP", and "TOB". Total observations available: 7.
#> Groups (n=4, named as 'order'):
#> [1] "Caryophanales"    "Enterobacterales" "Lactobacillales"  "Pseudomonadales" 
#> 


# }
```

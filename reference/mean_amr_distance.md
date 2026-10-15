# Calculate the Mean AMR Distance

Calculates a normalised mean for antimicrobial resistance between
multiple observations, to help to identify similar isolates without
comparing antibiograms by hand.

## Usage

``` r
mean_amr_distance(x, ...)

# S3 method for class 'sir'
mean_amr_distance(x, ..., combine_SI = TRUE)

# S3 method for class 'data.frame'
mean_amr_distance(x, ..., combine_SI = TRUE)

amr_distance_from_row(amr_distance, row)
```

## Arguments

- x:

  A vector of class [sir](https://amr-for-r.org/reference/as.sir.md),
  [mic](https://amr-for-r.org/reference/as.mic.md) or
  [disk](https://amr-for-r.org/reference/as.disk.md), or a
  [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns of any of these classes.

- ...:

  Variables to select. Supports [tidyselect
  language](https://tidyselect.r-lib.org/reference/starts_with.html)
  such as `where(is.mic)`, `starts_with(...)`, or `column1:column4`, and
  can thus also be [antimicrobial
  selectors](https://amr-for-r.org/reference/antimicrobial_selectors.md).

- combine_SI:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S, SDD, and I must be merged into one, so the input only
  consists of S+I vs. R (susceptible vs. resistant) - the default is
  `TRUE`.

- amr_distance:

  The outcome of `mean_amr_distance()`.

- row:

  An index, such as a row number.

## Details

The mean AMR distance is effectively [the
Z-score](https://en.wikipedia.org/wiki/Standard_score); a normalised
numeric value to compare AMR test results which can help to identify
similar isolates, without comparing antibiograms by hand.

MIC values (see [`as.mic()`](https://amr-for-r.org/reference/as.mic.md))
are transformed with [`log2()`](https://rdrr.io/r/base/Log.html) first;
their distance is thus calculated as
`(log2(x) - mean(log2(x))) / sd(log2(x))`.

SIR values (see [`as.sir()`](https://amr-for-r.org/reference/as.sir.md))
are transformed using `"S"` = 1, `"I"` = 2, and `"R"` = 3. If
`combine_SI` is `TRUE` (default), the `"I"` will be considered to be 1.

For data sets, the mean AMR distance will be calculated per column,
after which the mean per row will be returned, see *Examples*.

Use `amr_distance_from_row()` to subtract distances from the distance of
one row, see *Examples*.

## Interpretation

Isolates with distances less than 0.01 difference from each other should
be considered similar. Differences lower than 0.025 should be considered
suspicious.

## Examples

``` r
sir <- random_sir(10)
sir
#> Class <sir>
#>  [1] S I I I I I R R I R
mean_amr_distance(sir)
#>  [1] -0.621059 -0.621059 -0.621059 -0.621059 -0.621059 -0.621059  1.449138
#>  [8]  1.449138 -0.621059  1.449138

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] 0.5    0.008  0.001  0.0001 0.064  0.0002 0.064  0.0001 0.032  1     
mean_amr_distance(mic)
#>  [1]  1.198191624  0.000686808 -0.601499640 -1.268306334  0.602873256
#>  [6] -1.067577518  0.602873256 -1.268306334  0.402144440  1.398920440
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1]  1.198191624  0.000686808 -0.601499640 -1.268306334  0.602873256
#>  [6] -1.067577518  0.602873256 -1.268306334  0.402144440  1.398920440

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 38 37 42 30 40 24 45 36 49 44
mean_amr_distance(disk)
#>  [1] -0.06786705 -0.20360115  0.47506935 -1.15373986  0.20360115 -1.96814446
#>  [7]  0.88227166 -0.33933525  1.42520806  0.74653756

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    I   30   16  >=8
#> 2   B    R   30  <=8    4
#> 3   C    I   28  <=8    4
#> 4   D    I   24  <=8  >=8
#> 5   E    I   28  <=8    4
#> 6   F    S   30  <=8    4
#> 7   G    I   28  <=8    4
#> 8   H    S   31  <=8  >=8
#> 9   I    R   26  <=8    4
#> 10  J    I   28  <=8  >=8
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1]  1.08474449  0.40297922 -0.42682276 -0.41644967 -0.42682276 -0.18994784
#>  [7] -0.42682276  0.41261253 -0.07077061  0.05730016
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 2   B    R   30  <=8    4   -0.5454122
#> 3   C    I   28  <=8    4   -0.5454122
#> 5   E    I   28  <=8    4   -0.5454122
#> 6   F    S   30  <=8    4   -0.5454122
#> 7   G    I   28  <=8    4   -0.5454122
#> 9   I    R   26  <=8    4   -0.5454122
#> 4   D    I   24  <=8  >=8    0.4228336
#> 8   H    S   31  <=8  >=8    0.4228336
#> 10  J    I   28  <=8  >=8    0.4228336
#> 1   A    I   30   16  >=8    2.0039724

if (require("dplyr")) {
  y %>%
    mutate(
      amr_distance = mean_amr_distance(y),
      check_id_C = amr_distance_from_row(amr_distance, id == "C")
    ) %>%
    arrange(check_id_C)
}
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>    id amox cipr gent tobr amr_distance check_id_C
#> 1   C    I   28  <=8    4  -0.42682276 0.00000000
#> 2   E    I   28  <=8    4  -0.42682276 0.00000000
#> 3   G    I   28  <=8    4  -0.42682276 0.00000000
#> 4   D    I   24  <=8  >=8  -0.41644967 0.01037309
#> 5   F    S   30  <=8    4  -0.18994784 0.23687492
#> 6   I    R   26  <=8    4  -0.07077061 0.35605215
#> 7   J    I   28  <=8  >=8   0.05730016 0.48412292
#> 8   B    R   30  <=8    4   0.40297922 0.82980198
#> 9   H    S   31  <=8  >=8   0.41261253 0.83943529
#> 10  A    I   30   16  >=8   1.08474449 1.51156725
if (require("dplyr")) {
  # support for groups
  example_isolates %>%
    filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
    select(mo, TCY, carbapenems()) %>%
    group_by(mo) %>%
    mutate(dist = mean_amr_distance(.)) %>%
    arrange(mo, dist)
}
#> ℹ Using column mo as input for `mo_genus()`
#> ℹ Using column mo as input for `mo_species()`
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM (meropenem)
#> ℹ Calculating mean AMR distance based on columns "TCY", "IPM", and "MEM"
#> # A tibble: 63 × 5
#> # Groups:   mo [4]
#>    mo           TCY   IPM   MEM     dist
#>    <mo>         <sir> <sir> <sir>  <dbl>
#>  1 B_ENTRC_AVIM   S     S     NA   0    
#>  2 B_ENTRC_AVIM   S     S     NA   0    
#>  3 B_ENTRC_CSSL   NA    S     NA  NA    
#>  4 B_ENTRC_FACM   S     S     NA  -2.66 
#>  5 B_ENTRC_FACM   S     R     R   -0.423
#>  6 B_ENTRC_FACM   S     R     R   -0.423
#>  7 B_ENTRC_FACM   NA    R     R    0.224
#>  8 B_ENTRC_FACM   NA    R     R    0.224
#>  9 B_ENTRC_FACM   NA    R     R    0.224
#> 10 B_ENTRC_FACM   NA    R     R    0.224
#> # ℹ 53 more rows
```

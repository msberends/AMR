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
#> Class 'sir'
#>  [1] R R I R I I S R S I
mean_amr_distance(sir)
#>  [1]  1.1618950  1.1618950 -0.7745967  1.1618950 -0.7745967 -0.7745967
#>  [7] -0.7745967  1.1618950 -0.7745967 -0.7745967

mic <- random_mic(10)
mic
#> Class 'mic'
#>  [1] 0.032 0.5   1     >=8   4     0.016 >=8   1     0.004 0.008
mean_amr_distance(mic)
#>  [1] -0.7311752  0.2104422  0.4478776  1.1601837  0.9227483 -0.9686106
#>  [7]  1.1601837  0.4478776 -1.4434813 -1.2060459
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.7311752  0.2104422  0.4478776  1.1601837  0.9227483 -0.9686106
#>  [7]  1.1601837  0.4478776 -1.4434813 -1.2060459

disk <- random_disk(10)
disk
#> Class 'disk'
#>  [1] 50 49 38 33 31 17 42 43 46 37
mean_amr_distance(disk)
#>  [1]  1.15131286  1.05032051 -0.06059541 -0.56555720 -0.76754191 -2.18143490
#>  [7]  0.34337401  0.44436637  0.74734344 -0.16158777

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    I   27  >=2    8
#> 2   B    S   28    1    8
#> 3   C    R   33    1    8
#> 4   D    R   32    1   16
#> 5   E    I   25  0.5   16
#> 6   F    I   19  0.5    8
#> 7   G    S   23  0.5   16
#> 8   H    R   27  0.5    8
#> 9   I    S   29    1    8
#> 10  J    R   32  0.5   16
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent",
#>   and "tobr"
#>  [1]  0.08471751 -0.21572693  0.55391610  0.98093500 -0.26046456 -1.08721162
#>  [7] -0.37467261 -0.14625651 -0.15862291  0.62338653
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 6   F    I   19  0.5    8   -0.8163565
#> 8   H    R   27  0.5    8   -0.8163565
#> 2   B    S   28    1    8   -0.1012596
#> 3   C    R   33    1    8   -0.1012596
#> 9   I    S   29    1    8   -0.1012596
#> 5   E    I   25  0.5   16    0.1518893
#> 7   G    S   23  0.5   16    0.1518893
#> 10  J    R   32  0.5   16    0.1518893
#> 1   A    I   27  >=2    8    0.6138374
#> 4   D    R   32    1   16    0.8669863

if (require("dplyr")) {
  y %>%
    mutate(
      amr_distance = mean_amr_distance(y),
      check_id_C = amr_distance_from_row(amr_distance, id == "C")
    ) %>%
    arrange(check_id_C)
}
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent",
#>   and "tobr"
#>    id amox cipr gent tobr amr_distance check_id_C
#> 1   C    R   33    1    8   0.55391610 0.00000000
#> 2   J    R   32  0.5   16   0.62338653 0.06947042
#> 3   D    R   32    1   16   0.98093500 0.42701889
#> 4   A    I   27  >=2    8   0.08471751 0.46919859
#> 5   H    R   27  0.5    8  -0.14625651 0.70017262
#> 6   I    S   29    1    8  -0.15862291 0.71253901
#> 7   B    S   28    1    8  -0.21572693 0.76964304
#> 8   E    I   25  0.5   16  -0.26046456 0.81438066
#> 9   G    S   23  0.5   16  -0.37467261 0.92858871
#> 10  F    I   19  0.5    8  -1.08721162 1.64112773
if (require("dplyr")) {
  # support for groups
  example_isolates %>%
    filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
    select(mo, TCY, carbapenems()) %>%
    group_by(mo) %>%
    mutate(dist = mean_amr_distance(.)) %>%
    arrange(mo, dist)
}
#> ℹ Using column 'mo' as input for `mo_genus()`
#> ℹ Using column 'mo' as input for `mo_species()`
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
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

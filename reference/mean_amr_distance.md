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
#>  [1] S S I S I S I S I R
mean_amr_distance(sir)
#>  [1] -0.3162278 -0.3162278 -0.3162278 -0.3162278 -0.3162278 -0.3162278
#>  [7] -0.3162278 -0.3162278 -0.3162278  2.8460499

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] <=0.0001 0.0002   0.5      0.0005   0.004    0.001    0.016    0.0002  
#>  [9] 0.032    1       
mean_amr_distance(mic)
#>  [1] -1.13933320 -0.92906618  1.44436742 -0.65110831 -0.02030726 -0.44084129
#>  [7]  0.40022677 -0.92906618  0.61049379  1.65463443
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -1.13933320 -0.92906618  1.44436742 -0.65110831 -0.02030726 -0.44084129
#>  [7]  0.40022677 -0.92906618  0.61049379  1.65463443

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 49 30 42 50 16 27 12 13 28 32
mean_amr_distance(disk)
#>  [1]  1.37726856  0.00721083  0.87251045  1.44937686 -1.00230539 -0.20911407
#>  [7] -1.29073860 -1.21863029 -0.13700577  0.15142743

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    I   32    4    8
#> 2   B    R   32    4    8
#> 3   C    I   30    8    8
#> 4   D    S   16    4    8
#> 5   E    S   33    4   16
#> 6   F    R   32    4    8
#> 7   G    I   29    8    8
#> 8   H    I   31    8   16
#> 9   I    I   30    4   16
#> 10  J    R   20    8    8
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1] -0.351732344  0.165816825  0.045278388 -1.048629829  0.209372918
#>  [6]  0.165816825  0.001722296  0.606383651  0.078704640  0.127266630
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 1   A    I   32    4    8   -0.6978278
#> 2   B    R   32    4    8   -0.6978278
#> 4   D    S   16    4    8   -0.6978278
#> 6   F    R   32    4    8   -0.6978278
#> 3   C    I   30    8    8    0.2704180
#> 7   G    I   29    8    8    0.2704180
#> 10  J    R   20    8    8    0.2704180
#> 5   E    S   33    4   16    0.3372705
#> 9   I    I   30    4   16    0.3372705
#> 8   H    I   31    8   16    1.3055163

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
#> 1   C    I   30    8    8  0.045278388 0.00000000
#> 2   I    I   30    4   16  0.078704640 0.03342625
#> 3   G    I   29    8    8  0.001722296 0.04355609
#> 4   J    R   20    8    8  0.127266630 0.08198824
#> 5   B    R   32    4    8  0.165816825 0.12053844
#> 6   F    R   32    4    8  0.165816825 0.12053844
#> 7   E    S   33    4   16  0.209372918 0.16409453
#> 8   A    I   32    4    8 -0.351732344 0.39701073
#> 9   H    I   31    8   16  0.606383651 0.56110526
#> 10  D    S   16    4    8 -1.048629829 1.09390822
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

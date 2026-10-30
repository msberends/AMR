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
#>  [1] I R R I R I I R I I
mean_amr_distance(sir)
#>  [1] -0.7745967  1.1618950  1.1618950 -0.7745967  1.1618950 -0.7745967
#>  [7] -0.7745967  1.1618950 -0.7745967 -0.7745967

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] 0.0005 0.064  0.0005 0.064  0.0005 0.032  0.016  0.064  0.001  32    
mean_amr_distance(mic)
#>  [1] -0.98118944  0.42148922 -0.98118944  0.42148922 -0.98118944  0.22110656
#>  [7]  0.02072389  0.42148922 -0.78080678  2.21807699
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.98118944  0.42148922 -0.98118944  0.42148922 -0.98118944  0.22110656
#>  [7]  0.02072389  0.42148922 -0.78080678  2.21807699

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 50 24 48 45 42 45 39 47 39 18
mean_amr_distance(disk)
#>  [1]  0.9732102 -1.4834370  0.7842374  0.5007781  0.2173188  0.5007781
#>  [7] -0.0661405  0.6897510 -0.0661405 -2.0503556

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr   gent tobr
#> 1   A    S   24 <=0.25    8
#> 2   B    I   26 <=0.25    8
#> 3   C    S   21 <=0.25    8
#> 4   D    R   24 <=0.25    8
#> 5   E    I   26 <=0.25    8
#> 6   F    I   27    0.5    8
#> 7   G    S   21    0.5   16
#> 8   H    R   29 <=0.25   16
#> 9   I    R   26 <=0.25    8
#> 10  J    S   20 <=0.25    8
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1] -0.42631712 -0.25690939 -0.68042872  0.09123205 -0.25690939  0.42072154
#>  [7]  0.50542540  1.10767844  0.26063978 -0.76513259
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr   gent tobr amr_distance
#> 1   A    S   24 <=0.25    8   -0.4743416
#> 2   B    I   26 <=0.25    8   -0.4743416
#> 3   C    S   21 <=0.25    8   -0.4743416
#> 4   D    R   24 <=0.25    8   -0.4743416
#> 5   E    I   26 <=0.25    8   -0.4743416
#> 9   I    R   26 <=0.25    8   -0.4743416
#> 10  J    S   20 <=0.25    8   -0.4743416
#> 6   F    I   27    0.5    8    0.7115125
#> 8   H    R   29 <=0.25   16    0.7115125
#> 7   G    S   21    0.5   16    1.8973666

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
#>    id amox cipr   gent tobr amr_distance check_id_C
#> 1   C    S   21 <=0.25    8  -0.68042872 0.00000000
#> 2   J    S   20 <=0.25    8  -0.76513259 0.08470387
#> 3   A    S   24 <=0.25    8  -0.42631712 0.25411160
#> 4   B    I   26 <=0.25    8  -0.25690939 0.42351933
#> 5   E    I   26 <=0.25    8  -0.25690939 0.42351933
#> 6   D    R   24 <=0.25    8   0.09123205 0.77166077
#> 7   I    R   26 <=0.25    8   0.26063978 0.94106850
#> 8   F    I   27    0.5    8   0.42072154 1.10115026
#> 9   G    S   21    0.5   16   0.50542540 1.18585412
#> 10  H    R   29 <=0.25   16   1.10767844 1.78810716
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

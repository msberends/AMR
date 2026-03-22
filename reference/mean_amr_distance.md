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
#>  [1] S R S R I I R R R S
mean_amr_distance(sir)
#>  [1] -0.9486833  0.9486833 -0.9486833  0.9486833 -0.9486833 -0.9486833
#>  [7]  0.9486833  0.9486833  0.9486833 -0.9486833

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] 0.5      <=0.0001 0.25     0.0005   0.001    0.0002   8        <=0.0001
#>  [9] <=0.0001 0.004   
mean_amr_distance(mic)
#>  [1]  1.18999891 -0.86809178  1.02250714 -0.47918793 -0.31169616 -0.70060001
#>  [7]  1.85996599 -0.86809178 -0.86809178  0.02328738
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1]  1.18999891 -0.86809178  1.02250714 -0.47918793 -0.31169616 -0.70060001
#>  [7]  1.85996599 -0.86809178 -0.86809178  0.02328738

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 20 45 31 26 23 38 44 41 20 49
mean_amr_distance(disk)
#>  [1] -1.2414143  1.0239402 -0.2446583 -0.6977292 -0.9695717  0.3896410
#>  [7]  0.9333261  0.6614835 -1.2414143  1.3863970

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    R   33   16  <=1
#> 2   B    S   32  <=8   32
#> 3   C    S   27   16   32
#> 4   D    I   33  <=8    8
#> 5   E    I   22  <=8    2
#> 6   F    R   32  <=8    4
#> 7   G    S   31   16    8
#> 8   H    I   31  <=8  <=1
#> 9   I    R   28   16    2
#> 10  J    S   22   16    4
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1]  0.52677313  0.16501937  0.34372779 -0.05155946 -0.97765805  0.26901032
#>  [7]  0.30452889 -0.58337098  0.36899264 -0.36546366
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 8   H    I   31  <=8  <=1   -1.0808937
#> 5   E    I   22  <=8    2   -0.8051882
#> 6   F    R   32  <=8    4   -0.5294827
#> 4   D    I   33  <=8    8   -0.2537773
#> 1   A    R   33   16  <=1   -0.1322104
#> 9   I    R   28   16    2    0.1434951
#> 2   B    S   32  <=8   32    0.2976337
#> 10  J    S   22   16    4    0.4192006
#> 7   G    S   31   16    8    0.6949060
#> 3   C    S   27   16   32    1.2463170

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
#> 1   C    S   27   16   32   0.34372779 0.00000000
#> 2   I    R   28   16    2   0.36899264 0.02526485
#> 3   G    S   31   16    8   0.30452889 0.03919890
#> 4   F    R   32  <=8    4   0.26901032 0.07471748
#> 5   B    S   32  <=8   32   0.16501937 0.17870842
#> 6   A    R   33   16  <=1   0.52677313 0.18304534
#> 7   D    I   33  <=8    8  -0.05155946 0.39528726
#> 8   J    S   22   16    4  -0.36546366 0.70919145
#> 9   H    I   31  <=8  <=1  -0.58337098 0.92709877
#> 10  E    I   22  <=8    2  -0.97765805 1.32138584
if (require("dplyr")) {
  # support for groups
  example_isolates %>%
    filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
    select(mo, TCY, carbapenems()) %>%
    group_by(mo) %>%
    mutate(dist = mean_amr_distance(.)) %>%
    arrange(mo, dist)
}
#> ℹ Using column 'mo' as input for `?mo_genus()`
#> ℹ Using column 'mo' as input for `?mo_species()`
#> ℹ For `?carbapenems()` using columns IPM (imipenem) and MEM (meropenem)
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

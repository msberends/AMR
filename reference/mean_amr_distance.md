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
#>  [1] R S I R R I S R I R
mean_amr_distance(sir)
#>  [1]  0.9486833 -0.9486833 -0.9486833  0.9486833  0.9486833 -0.9486833
#>  [7] -0.9486833  0.9486833 -0.9486833  0.9486833

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] 0.0002   0.032    0.002    8        <=0.0001 0.001    0.5      <=0.0001
#>  [9] 0.002    0.002   
mean_amr_distance(mic)
#>  [1] -0.8220689  0.5401854 -0.2040199  2.0222300 -1.0081203 -0.3900712
#>  [7]  1.2780248 -1.0081203 -0.2040199 -0.2040199
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.8220689  0.5401854 -0.2040199  2.0222300 -1.0081203 -0.3900712
#>  [7]  1.2780248 -1.0081203 -0.2040199 -0.2040199

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 49 26 48 42 31 50 13 43 50 35
mean_amr_distance(disk)
#>  [1]  0.8341912 -1.0285658  0.7532018  0.2672651 -0.6236187  0.9151806
#>  [7] -2.0814285  0.3482546  0.9151806 -0.2996609

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    S   24    4    2
#> 2   B    I   27    4    2
#> 3   C    S   18    4    2
#> 4   D    S   30    4  <=1
#> 5   E    I   23    4    2
#> 6   F    I   19    4    2
#> 7   G    R   23    4  <=1
#> 8   H    I   28    4  <=1
#> 9   I    R   28    4  <=1
#> 10  J    I   26    4  <=1
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1]  0.10747691  0.36066176 -0.39889277 -0.01860893  0.02308197 -0.31449782
#>  [7]  0.18119585 -0.18739883  0.60317059 -0.35618872
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 4   D    S   30    4  <=1   -0.9486833
#> 7   G    R   23    4  <=1   -0.9486833
#> 8   H    I   28    4  <=1   -0.9486833
#> 9   I    R   28    4  <=1   -0.9486833
#> 10  J    I   26    4  <=1   -0.9486833
#> 1   A    S   24    4    2    0.9486833
#> 2   B    I   27    4    2    0.9486833
#> 3   C    S   18    4    2    0.9486833
#> 5   E    I   23    4    2    0.9486833
#> 6   F    I   19    4    2    0.9486833

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
#> 1   C    S   18    4    2  -0.39889277 0.00000000
#> 2   J    I   26    4  <=1  -0.35618872 0.04270405
#> 3   F    I   19    4    2  -0.31449782 0.08439495
#> 4   H    I   28    4  <=1  -0.18739883 0.21149394
#> 5   D    S   30    4  <=1  -0.01860893 0.38028384
#> 6   E    I   23    4    2   0.02308197 0.42197474
#> 7   A    S   24    4    2   0.10747691 0.50636968
#> 8   G    R   23    4  <=1   0.18119585 0.58008862
#> 9   B    I   27    4    2   0.36066176 0.75955453
#> 10  I    R   28    4  <=1   0.60317059 1.00206336
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

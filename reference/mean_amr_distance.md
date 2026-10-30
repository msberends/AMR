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
#>  [1] I I S I R I S S S R
mean_amr_distance(sir)
#>  [1] -0.4743416 -0.4743416 -0.4743416 -0.4743416  1.8973666 -0.4743416
#>  [7] -0.4743416 -0.4743416 -0.4743416  1.8973666

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] <=0.0001 64       0.001    64       0.0002   <=0.0001 0.008    0.064   
#>  [9] 0.0002   0.125   
mean_amr_distance(mic)
#>  [1] -0.9348351  1.6473085 -0.4901118  1.6473085 -0.8009600 -0.9348351
#>  [7] -0.0884866  0.3131386 -0.8009600  0.4424330
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.9348351  1.6473085 -0.4901118  1.6473085 -0.8009600 -0.9348351
#>  [7] -0.0884866  0.3131386 -0.8009600  0.4424330

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 43 46 40 43 47 47 50 20 46 31
mean_amr_distance(disk)
#>  [1]  0.1859655  0.5141399 -0.1422089  0.1859655  0.6235314  0.6235314
#>  [7]  0.9517058 -2.3300383  0.5141399 -1.1267321

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    I   31   16    4
#> 2   B    S   30  <=8  <=2
#> 3   C    S   27  <=8  <=2
#> 4   D    S   32  <=8  <=2
#> 5   E    R   27  <=8  <=2
#> 6   F    I   32  <=8    4
#> 7   G    R   26  <=8  <=2
#> 8   H    R   30   16    4
#> 9   I    I   19  <=8  <=2
#> 10  J    S   19  <=8  <=2
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1]  0.87190417 -0.29007013 -0.44456433 -0.18707400  0.07298484  0.33047517
#>  [7]  0.02148677  1.33795527 -0.85654887 -0.85654887
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 2   B    S   30  <=8  <=2   -0.5477003
#> 3   C    S   27  <=8  <=2   -0.5477003
#> 4   D    S   32  <=8  <=2   -0.5477003
#> 5   E    R   27  <=8  <=2   -0.5477003
#> 7   G    R   26  <=8  <=2   -0.5477003
#> 9   I    I   19  <=8  <=2   -0.5477003
#> 10  J    S   19  <=8  <=2   -0.5477003
#> 6   F    I   32  <=8    4    0.4873980
#> 1   A    I   31   16    4    1.6732521
#> 8   H    R   30   16    4    1.6732521

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
#> 1   C    S   27  <=8  <=2  -0.44456433  0.0000000
#> 2   B    S   30  <=8  <=2  -0.29007013  0.1544942
#> 3   D    S   32  <=8  <=2  -0.18707400  0.2574903
#> 4   I    I   19  <=8  <=2  -0.85654887  0.4119845
#> 5   J    S   19  <=8  <=2  -0.85654887  0.4119845
#> 6   G    R   26  <=8  <=2   0.02148677  0.4660511
#> 7   E    R   27  <=8  <=2   0.07298484  0.5175492
#> 8   F    I   32  <=8    4   0.33047517  0.7750395
#> 9   A    I   31   16    4   0.87190417  1.3164685
#> 10  H    R   30   16    4   1.33795527  1.7825196
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

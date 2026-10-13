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
#>  [1] S R I R R R R I I I
mean_amr_distance(sir)
#>  [1] -0.9486833  0.9486833 -0.9486833  0.9486833  0.9486833  0.9486833
#>  [7]  0.9486833 -0.9486833 -0.9486833 -0.9486833

mic <- random_mic(10)
mic
#> Class <mic>
#>  [1] 0.032    0.032    0.125    0.25     0.016    0.004    <=0.0005 0.016   
#>  [9] 0.25     0.5     
mean_amr_distance(mic)
#>  [1] -0.02823448 -0.02823448  0.61480238  0.94191703 -0.35534914 -1.00957846
#>  [7] -1.99092244 -0.35534914  0.94191703  1.26903169
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.02823448 -0.02823448  0.61480238  0.94191703 -0.35534914 -1.00957846
#>  [7] -1.99092244 -0.35534914  0.94191703  1.26903169

disk <- random_disk(10)
disk
#> Class <disk>
#>  [1] 48 45 48 40 44  9 39 49 39 29
mean_amr_distance(disk)
#>  [1]  0.74202711  0.49468474  0.74202711  0.08244746  0.41223728 -2.47342369
#>  [7]  0.00000000  0.82447456  0.00000000 -0.82447456

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    R   26  <=2  >=2
#> 2   B    R   28  <=2    1
#> 3   C    R   30  <=2    1
#> 4   D    R   30  <=2    1
#> 5   E    I   31    4    1
#> 6   F    S   28    4    1
#> 7   G    S   32  <=2  >=2
#> 8   H    I   25  <=2    1
#> 9   I    S   21  <=2    1
#> 10  J    R   20    4    1
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent", and
#>   "tobr"
#>  [1]  0.48908793  0.01826958  0.14037829  0.14037829  0.24464016  0.06147710
#>  [7]  0.38107241 -0.63923513 -0.88345255  0.04738391
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 2   B    R   28  <=2    1   -0.5477003
#> 3   C    R   30  <=2    1   -0.5477003
#> 4   D    R   30  <=2    1   -0.5477003
#> 8   H    I   25  <=2    1   -0.5477003
#> 9   I    S   21  <=2    1   -0.5477003
#> 5   E    I   31    4    1    0.4873980
#> 6   F    S   28    4    1    0.4873980
#> 10  J    R   20    4    1    0.4873980
#> 1   A    R   26  <=2  >=2    0.6381538
#> 7   G    S   32  <=2  >=2    0.6381538

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
#> 1   C    R   30  <=2    1   0.14037829 0.00000000
#> 2   D    R   30  <=2    1   0.14037829 0.00000000
#> 3   F    S   28    4    1   0.06147710 0.07890119
#> 4   J    R   20    4    1   0.04738391 0.09299437
#> 5   E    I   31    4    1   0.24464016 0.10426187
#> 6   B    R   28  <=2    1   0.01826958 0.12210871
#> 7   G    S   32  <=2  >=2   0.38107241 0.24069412
#> 8   A    R   26  <=2  >=2   0.48908793 0.34870964
#> 9   H    I   25  <=2    1  -0.63923513 0.77961342
#> 10  I    S   21  <=2    1  -0.88345255 1.02383084
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

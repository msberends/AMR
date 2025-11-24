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
#>  [1] I I R I R S S S I S
mean_amr_distance(sir)
#>  [1] -0.4743416 -0.4743416  1.8973666 -0.4743416  1.8973666 -0.4743416
#>  [7] -0.4743416 -0.4743416 -0.4743416 -0.4743416

mic <- random_mic(10)
mic
#> Class 'mic'
#>  [1] 0.004  2      0.002  0.0001 0.004  0.002  >=4    0.0002 0.032  0.004 
mean_amr_distance(mic)
#>  [1] -0.2047915  1.5799751 -0.4038557 -1.2641969 -0.2047915 -0.4038557
#>  [7]  1.7790393 -1.0651327  0.3924011 -0.2047915
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1] -0.2047915  1.5799751 -0.4038557 -1.2641969 -0.2047915 -0.4038557
#>  [7]  1.7790393 -1.0651327  0.3924011 -0.2047915

disk <- random_disk(10)
disk
#> Class 'disk'
#>  [1] 43 12 28 32 22 31 35 25 43 35
mean_amr_distance(disk)
#>  [1]  1.30998909 -1.96498364 -0.27467513  0.14790199 -0.90854082  0.04225771
#>  [7]  0.46483484 -0.59160798  1.30998909  0.46483484

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    S   31    2 >=16
#> 2   B    S   27  <=1    8
#> 3   C    R   25    2    4
#> 4   D    R   25  <=1    2
#> 5   E    I   31  <=1    2
#> 6   F    S   32  <=1    8
#> 7   G    I   29    2    2
#> 8   H    S   18  <=1    4
#> 9   I    S   28  <=1    4
#> 10  J    R   17  <=1    2
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent",
#>   and "tobr"
#>  [1]  0.90606144 -0.03989270  0.66241774 -0.09230226 -0.32300020  0.19914999
#>  [7]  0.09893189 -0.70734036 -0.22925499 -0.47477055
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 4   D    R   25  <=1    2   -0.7848712
#> 5   E    I   31  <=1    2   -0.7848712
#> 10  J    R   17  <=1    2   -0.7848712
#> 8   H    S   18  <=1    4   -0.3105295
#> 9   I    S   28  <=1    4   -0.3105295
#> 2   B    S   27  <=1    8    0.1638121
#> 6   F    S   32  <=1    8    0.1638121
#> 7   G    I   29    2    2    0.2502272
#> 3   C    R   25    2    4    0.7245688
#> 1   A    S   31    2 >=16    1.6732521

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
#> 1   C    R   25    2    4   0.66241774  0.0000000
#> 2   A    S   31    2 >=16   0.90606144  0.2436437
#> 3   F    S   32  <=1    8   0.19914999  0.4632678
#> 4   G    I   29    2    2   0.09893189  0.5634858
#> 5   B    S   27  <=1    8  -0.03989270  0.7023104
#> 6   D    R   25  <=1    2  -0.09230226  0.7547200
#> 7   I    S   28  <=1    4  -0.22925499  0.8916727
#> 8   E    I   31  <=1    2  -0.32300020  0.9854179
#> 9   J    R   17  <=1    2  -0.47477055  1.1371883
#> 10  H    S   18  <=1    4  -0.70734036  1.3697581
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

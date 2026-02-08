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
#>  [1] R I R I I S R S I I
mean_amr_distance(sir)
#>  [1]  1.449138 -0.621059  1.449138 -0.621059 -0.621059 -0.621059  1.449138
#>  [8] -0.621059 -0.621059 -0.621059

mic <- random_mic(10)
mic
#> Class 'mic'
#>  [1] 0.25     0.5      8        2        0.004    4        1        0.001   
#>  [9] 0.004    <=0.0005
mean_amr_distance(mic)
#>  [1]  0.2626139  0.4520471  1.2097801  0.8309136 -0.8675039  1.0203468
#>  [7]  0.6414804 -1.2463704 -0.8675039 -1.4358036
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))
#>  [1]  0.2626139  0.4520471  1.2097801  0.8309136 -0.8675039  1.0203468
#>  [7]  0.6414804 -1.2463704 -0.8675039 -1.4358036

disk <- random_disk(10)
disk
#> Class 'disk'
#>  [1] 49 38 33 31 17 42 43 46 37 46
mean_amr_distance(disk)
#>  [1]  1.14152462 -0.02113934 -0.54962296 -0.76101641 -2.24077054  0.40164755
#>  [7]  0.50734427  0.82443445 -0.12683607  0.82443445

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
#>    id amox cipr gent tobr
#> 1   A    S   28  >=2    8
#> 2   B    R   33    1    8
#> 3   C    R   32    1    8
#> 4   D    I   25    1   16
#> 5   E    I   19  0.5   16
#> 6   F    S   23  0.5    8
#> 7   G    R   27  0.5   16
#> 8   H    S   29  0.5    8
#> 9   I    R   32    1    8
#> 10  J    R   32  0.5   16
mean_amr_distance(y)
#> ℹ Calculating mean AMR distance based on columns "amox", "cipr", "gent",
#>   and "tobr"
#>  [1]  0.06974787  0.45859464  0.40418392  0.03309016 -0.65092262 -0.91740267
#>  [7]  0.25870477 -0.59093836  0.40418392  0.53075837
y$amr_distance <- mean_amr_distance(y, is.mic(y))
#> ℹ Calculating mean AMR distance based on columns "gent" and "tobr"
y[order(y$amr_distance), ]
#>    id amox cipr gent tobr amr_distance
#> 6   F    S   23  0.5    8   -0.8163565
#> 8   H    S   29  0.5    8   -0.8163565
#> 2   B    R   33    1    8   -0.1012596
#> 3   C    R   32    1    8   -0.1012596
#> 9   I    R   32    1    8   -0.1012596
#> 5   E    I   19  0.5   16    0.1518893
#> 7   G    R   27  0.5   16    0.1518893
#> 10  J    R   32  0.5   16    0.1518893
#> 1   A    S   28  >=2    8    0.6138374
#> 4   D    I   25    1   16    0.8669863

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
#> 1   C    R   32    1    8   0.40418392 0.00000000
#> 2   I    R   32    1    8   0.40418392 0.00000000
#> 3   B    R   33    1    8   0.45859464 0.05441072
#> 4   J    R   32  0.5   16   0.53075837 0.12657445
#> 5   G    R   27  0.5   16   0.25870477 0.14547915
#> 6   A    S   28  >=2    8   0.06974787 0.33443605
#> 7   D    I   25    1   16   0.03309016 0.37109376
#> 8   H    S   29  0.5    8  -0.59093836 0.99512228
#> 9   E    I   19  0.5   16  -0.65092262 1.05510655
#> 10  F    S   23  0.5    8  -0.91740267 1.32158659
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

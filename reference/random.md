# Random MIC Values/Disk Zones/SIR Generation

These functions can be used for generating random MIC values and disk
diffusion diameters, for AMR data analysis practice. By providing a
microorganism and antimicrobial drug, the generated results will reflect
reality as much as possible.

## Usage

``` r
random_mic(
  size = NULL,
  mo = NULL,
  ab = NULL,
  skew = "right",
  severity = 1,
  ...
)

random_disk(
  size = NULL,
  mo = NULL,
  ab = NULL,
  skew = "left",
  severity = 1,
  ...
)

random_sir(size = NULL, prob_SIR = c(0.33, 0.33, 0.33), ...)
```

## Arguments

- size:

  Desired size of the returned vector. If used in a
  [data.frame](https://rdrr.io/r/base/data.frame.html) call or `dplyr`
  verb, will get the current (group) size if left blank.

- mo:

  Any [character](https://rdrr.io/r/base/character.html) that can be
  coerced to a valid microorganism code with
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md). Can be the same
  length as `size`.

- ab:

  Any [character](https://rdrr.io/r/base/character.html) that can be
  coerced to a valid antimicrobial drug code with
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

- skew:

  Direction of skew for MIC or disk values, either `"right"` or
  `"left"`. A left-skewed distribution has the majority of the data on
  the right.

- severity:

  Skew severity; higher values will increase the skewedness. Default is
  `2`; use `0` to prevent skewedness.

- ...:

  Ignored, only in place to allow future extensions.

- prob_SIR:

  A vector of length 3: the probabilities for "S" (1st value), "I" (2nd
  value) and "R" (3rd value).

## Value

class `mic` for `random_mic()` (see
[`as.mic()`](https://amr-for-r.org/reference/as.mic.md)) and class
`disk` for `random_disk()` (see
[`as.disk()`](https://amr-for-r.org/reference/as.disk.md))

## Details

Internally, MIC and disk zone values are sampled based on clinical
breakpoints defined in the
[clinical_breakpoints](https://amr-for-r.org/reference/clinical_breakpoints.md)
data set. To create specific generated values per bug or drug, set the
`mo` and/or `ab` argument. The MICs are sampled on a log2 scale and
disks linearly, using weighted probabilities. The weights are based on
the `skew` and `severity` arguments:

- `skew = "right"` places more emphasis on lower MIC or higher disk
  values.

- `skew = "left"` places more emphasis on higher MIC or lower disk
  values.

- `severity` controls the exponential bias applied.

## Examples

``` r
random_mic(25)
#> Class <mic>
#>  [1] 0.25   0.064  0.5    0.004  0.0001 0.032  0.008  16     0.016  0.125 
#> [11] 0.0002 0.001  0.008  0.004  0.064  256    0.5    32     512    0.0002
#> [21] 0.004  1      0.002  0.0005 0.064 
random_disk(25)
#> Class <disk>
#>  [1] 49 44 31 11 27 33 39 24 19 43 30 44 33 31 22 13 37 19 26 29 44 12 46 18 27
random_sir(25)
#> Class <sir>
#>  [1] I R I S R S R R I R R I S S S I R R S I R S R I R

# add more skewedness, make more realistic by setting a bug and/or drug:
disks <- random_disk(100, severity = 2, mo = "Escherichia coli", ab = "CIP")
plot(disks)

# `plot()` and `ggplot2::autoplot()` allow for coloured bars if `mo` and `ab` are set
plot(disks, mo = "Escherichia coli", ab = "CIP", guideline = "CLSI 2025")


# \donttest{
random_mic(25, "Klebsiella pneumoniae") # range 0.0625-64
#> Class <mic>
#>  [1] 0.004  0.0002 0.008  0.25   0.0005 0.0002 8      0.032  0.0005 0.016 
#> [11] 0.032  32     1      0.0005 0.016  0.002  0.016  0.0005 0.004  0.25  
#> [21] 0.0005 0.016  0.001  0.25   0.032 
random_mic(25, "Klebsiella pneumoniae", "meropenem") # range 0.0625-16
#> Class <mic>
#>  [1] 0.5 >=1 0.5 0.5 0.5 0.5 >=1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 >=1 0.5 >=1 0.5 0.5
#> [20] 0.5 0.5 >=1 >=1 >=1 >=1
random_mic(25, "Streptococcus pneumoniae", "meropenem") # range 0.0625-4
#> Class <mic>
#>  [1] 0.25    <=0.125 <=0.125 0.25    0.25    <=0.125 0.25    0.25    0.25   
#> [10] 0.25    <=0.125 <=0.125 0.25    0.25    <=0.125 <=0.125 <=0.125 0.25   
#> [19] <=0.125 <=0.125 <=0.125 <=0.125 0.25    0.25    0.25   

random_disk(25, "Klebsiella pneumoniae") # range 8-50
#> Class <disk>
#>  [1] 10 32 34 18 17 31 31 21 18 19 27 18 20 11 34 30 28 20 29  9 21 33 34 25 17
random_disk(25, "Klebsiella pneumoniae", "ampicillin") # range 11-17
#> Class <disk>
#>  [1] 17 17 20 22 13 21 21 14 11 18 18 18 21 18 14 22 19 20 18 16 17 12 22 21 19
random_disk(25, "Streptococcus pneumoniae", "ampicillin") # range 12-27
#> Class <disk>
#>  [1] 24 31 26 26 33 16 29 17 28 27 27 35 33 33 29 32 18 33 33 29 28 32 25 17 31
# }
```

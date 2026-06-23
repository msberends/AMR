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
#>  [1] 1      0.032  0.064  1      0.25   0.125  >=128  0.0002 32     0.25  
#> [11] 0.008  0.125  0.001  0.032  8      >=128  0.032  0.008  8      >=128 
#> [21] 0.125  0.0001 0.5    64     1     
random_disk(25)
#> Class <disk>
#>  [1] 49 39 29 30 16 23 44 50 25 28 29 14 44 34  8 33 33 34 38 42 10 42 43 11 46
random_sir(25)
#> Class <sir>
#>  [1] I S S S R R S R I I I I S R S I S S S R R I R S I

# add more skewedness, make more realistic by setting a bug and/or drug:
disks <- random_disk(100, severity = 2, mo = "Escherichia coli", ab = "CIP")
plot(disks)

# `plot()` and `ggplot2::autoplot()` allow for coloured bars if `mo` and `ab` are set
plot(disks, mo = "Escherichia coli", ab = "CIP", guideline = "CLSI 2025")


# \donttest{
random_mic(25, "Klebsiella pneumoniae") # range 0.0625-64
#> Class <mic>
#>  [1] 256      <=0.0001 <=0.0001 0.016    1        0.016    2        8       
#>  [9] 0.5      0.004    0.032    0.004    0.001    0.001    0.064    0.002   
#> [17] 0.032    0.004    0.008    16       0.125    <=0.0001 0.001    0.008   
#> [25] 0.001   
random_mic(25, "Klebsiella pneumoniae", "meropenem") # range 0.0625-16
#> Class <mic>
#>  [1] 0.5    <=0.25 0.5    4      1      1      0.5    0.5    0.5    <=0.25
#> [11] 0.5    0.5    2      1      <=0.25 <=0.25 8      <=0.25 1      1     
#> [21] <=0.25 0.5    16     2      0.5   
random_mic(25, "Streptococcus pneumoniae", "meropenem") # range 0.0625-4
#> Class <mic>
#>  [1] 0.125  0.125  0.125  0.125  0.125  >=0.25 0.125  0.125  0.125  0.125 
#> [11] >=0.25 0.125  >=0.25 0.125  0.125  0.125  >=0.25 0.125  >=0.25 >=0.25
#> [21] >=0.25 0.125  0.125  0.125  >=0.25

random_disk(25, "Klebsiella pneumoniae") # range 8-50
#> Class <disk>
#>  [1] 18 32 21 19 23 10 30 18 22 34 19 32 28 18 25  7 27 23 28 17 29 18 22 28 23
random_disk(25, "Klebsiella pneumoniae", "ampicillin") # range 11-17
#> Class <disk>
#>  [1] 22 18 13 22 17 12 20 20 16 20 22 22 18 20 20 22 22 16 13 17 15 17 16 22 15
random_disk(25, "Streptococcus pneumoniae", "ampicillin") # range 12-27
#> Class <disk>
#>  [1] 31 30 30 28 24 35 28 34 28 29 27 29 24 21 13 16 34 19 26 20 27 32 32 24 26
# }
```

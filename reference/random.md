# Generate Random MIC Values, Disk Zones, or SIR Results

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
  guideline = getOption("AMR_guideline", "EUCAST"),
  ...
)

random_disk(
  size = NULL,
  mo = NULL,
  ab = NULL,
  skew = "left",
  severity = 1,
  guideline = getOption("AMR_guideline", "EUCAST"),
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

- guideline:

  A guideline name to base the generated values on, either `"EUCAST"`
  (default) or `"CLSI"` for the latest implemented version, or a
  guideline including a year (e.g., `"CLSI 2019"`). The default can be
  changed with the package option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), see
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md).

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

Internally, MIC and disk zone values are sampled based on the human
clinical breakpoints of `guideline` in the
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
#>  [1] 0.032    0.002    0.25     0.064    0.064    0.125    0.004    4       
#>  [9] 0.008    2        0.001    0.016    0.001    128      <=0.0005 512     
#> [17] 0.001    256      4        4        0.004    0.002    0.125    0.004   
#> [25] 0.001   
random_disk(25)
#> Class <disk>
#>  [1] 29 36 36 29 43 10 39 39 39 49 21 43 34 45 38 43 33 47 48 26 26 39 27 29 26
random_sir(25)
#> Class <sir>
#>  [1] I I R S R S S R I S S R S R S I S S I S R I R I R

# add more skewedness, make more realistic by setting a bug and/or drug:
disks <- random_disk(100, severity = 2, mo = "Escherichia coli", ab = "CIP")
plot(disks)

# `plot()` and `ggplot2::autoplot()` allow for coloured bars if `mo` and `ab` are set
plot(disks, mo = "Escherichia coli", ab = "CIP", guideline = "CLSI")


# \donttest{
random_mic(25, "Klebsiella pneumoniae") # range 0.0625-64
#> Class <mic>
#>  [1] 0.064    0.008    0.25     0.0005   0.016    0.25     0.032    0.0005  
#>  [9] 0.016    1        0.002    0.0002   0.125    8        0.001    0.032   
#> [17] 2        2        0.0005   <=0.0001 0.125    <=0.0001 0.0005   0.004   
#> [25] 0.002   
random_mic(25, "Klebsiella pneumoniae", "meropenem") # range 0.0625-16
#> Class <mic>
#>  [1] <=2 8   8   8   8   <=2 <=2 8   <=2 <=2 <=2 8   <=2 <=2 4   16  4   4   16 
#> [20] <=2 16  8   4   <=2 4  
random_mic(25, "Streptococcus pneumoniae", "meropenem") # range 0.0625-4
#> Class <mic>
#>  [1] 0.5     0.5     0.25    <=0.064 0.25    1       <=0.064 <=0.064 2      
#> [10] 0.25    0.5     1       0.125   <=0.064 0.125   0.25    0.5     <=0.064
#> [19] 1       <=0.064 <=0.064 2       <=0.064 0.125   0.125  

random_disk(25, "Klebsiella pneumoniae") # range 8-50
#> Class <disk>
#>  [1] 31 28 32 34 27 31 24 28 29 27 23 33 29 27 26 32 34 33 34 19 23 30 34 21 33
random_disk(25, "Klebsiella pneumoniae", "ampicillin") # range 11-17
#> Class <disk>
#>  [1] 18 15 13 17 17 13 18 17 18 17 16 18 18 14 18 18 13 14 18 16 15 18 18 14 11
random_disk(25, "Streptococcus pneumoniae", "ampicillin") # range 12-27
#> Class <disk>
#>  [1] 22 28 21 17 25 25 15 25 22 27 20 25 24 19 22 22 28 28 24 17 22 28 26 25 28
# }
```

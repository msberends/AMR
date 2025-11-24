# Transform Input to Disk Diffusion Diameters

This transforms a vector to a new class `disk`, which is a disk
diffusion growth zone size (around an antibiotic disk) in millimetres
between 0 and 50.

## Usage

``` r
as.disk(x, na.rm = FALSE)

NA_disk_

is.disk(x)
```

## Format

An object of class `disk` (inherits from `integer`) of length 1.

## Arguments

- x:

  Vector.

- na.rm:

  A [logical](https://rdrr.io/r/base/logical.html) indicating whether
  missing values should be removed.

## Value

An [integer](https://rdrr.io/r/base/integer.html) with additional class
`disk`

## Details

Interpret disk values as SIR values with
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md). It supports
guidelines from EUCAST and CLSI.

Disk diffusion growth zone sizes must be between 0 and 50 millimetres.
Values higher than 50 but lower than 100 will be maximised to 50. All
others input values outside the 0-50 range will return `NA`.

`NA_disk_` is a missing value of the new `disk` class.

## See also

[`as.sir()`](https://amr-for-r.org/reference/as.sir.md)

## Examples

``` r
# transform existing disk zones to the `disk` class (using base R)
df <- data.frame(
  microorganism = "Escherichia coli",
  AMP = 20,
  CIP = 14,
  GEN = 18,
  TOB = 16
)
df[, 2:5] <- lapply(df[, 2:5], as.disk)
str(df)
#> 'data.frame':    1 obs. of  5 variables:
#>  $ microorganism: chr "Escherichia coli"
#>  $ AMP          : 'disk' int 20
#>  $ CIP          : 'disk' int 14
#>  $ GEN          : 'disk' int 18
#>  $ TOB          : 'disk' int 16

# \donttest{
# transforming is easier with dplyr:
if (require("dplyr")) {
  df %>% mutate(across(AMP:TOB, as.disk))
}
#>      microorganism AMP CIP GEN TOB
#> 1 Escherichia coli  20  14  18  16
# }

# interpret disk values, see ?as.sir
as.sir(
  x = as.disk(18),
  mo = "Strep pneu", # `mo` will be coerced with as.mo()
  ab = "ampicillin", # and `ab` with as.ab()
  guideline = "EUCAST"
)
#> Class 'sir'
#> [1] R

# interpret whole data set, pretend to be all from urinary tract infections:
as.sir(df, uti = TRUE)
#>      microorganism AMP  CIP GEN TOB
#> 1 Escherichia coli   S <NA>   S   S
```

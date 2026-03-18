# Age in Years of Individuals

Calculates age in years based on a reference date, which is the system
date at default.

## Usage

``` r
age(x, reference = Sys.Date(), exact = FALSE, na.rm = FALSE, ...)
```

## Arguments

- x:

  Date(s), [character](https://rdrr.io/r/base/character.html) (vectors)
  will be coerced with
  [`as.POSIXlt()`](https://rdrr.io/r/base/as.POSIXlt.html).

- reference:

  Reference date(s) (default is today),
  [character](https://rdrr.io/r/base/character.html) (vectors) will be
  coerced with [`as.POSIXlt()`](https://rdrr.io/r/base/as.POSIXlt.html).

- exact:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  age calculation should be exact, i.e. with decimals. It divides the
  number of days of
  [year-to-date](https://en.wikipedia.org/wiki/Year-to-date) (YTD) of
  `x` by the number of days in the year of `reference` (either 365 or
  366).

- na.rm:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  missing values should be removed.

- ...:

  Arguments passed on to
  [`as.POSIXlt()`](https://rdrr.io/r/base/as.POSIXlt.html), such as
  `origin`.

## Value

An [integer](https://rdrr.io/r/base/integer.html) (no decimals) if
`exact = FALSE`, a [double](https://rdrr.io/r/base/double.html) (with
decimals) otherwise

## Details

Ages below 0 will be returned as `NA` with a warning. Ages above 120
will only give a warning.

This function vectorises over both `x` and `reference`, meaning that
either can have a length of 1 while the other argument has a larger
length.

## See also

To split ages into groups, use the
[`age_groups()`](https://amr-for-r.org/reference/age_groups.md)
function.

## Examples

``` r
# 10 random pre-Y2K birth dates
df <- data.frame(birth_date = as.Date("2000-01-01") - runif(10) * 25000)

# add ages
df$age <- age(df$birth_date)

# add exact ages
df$age_exact <- age(df$birth_date, exact = TRUE)

# add age at millenium switch
df$age_at_y2k <- age(df$birth_date, "2000-01-01")

df
#>    birth_date age age_exact age_at_y2k
#> 1  1999-06-30  26  26.71507          0
#> 2  1968-01-29  58  58.13151         31
#> 3  1965-12-05  60  60.28219         34
#> 4  1980-03-01  46  46.04658         19
#> 5  1949-11-01  76  76.37534         50
#> 6  1947-02-14  79  79.08767         52
#> 7  1940-02-19  86  86.07397         59
#> 8  1988-01-10  38  38.18356         11
#> 9  1997-08-27  28  28.55616          2
#> 10 1978-01-26  48  48.13973         21
```

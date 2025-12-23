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
#> 1  1980-02-27  45  45.81918         19
#> 2  1953-07-26  72  72.41096         46
#> 3  1949-09-02  76  76.30685         50
#> 4  1986-08-03  39  39.38904         13
#> 5  1932-11-19  93  93.09315         67
#> 6  1949-03-30  76  76.73425         50
#> 7  1996-06-23  29  29.50137          3
#> 8  1963-09-16  62  62.26849         36
#> 9  1952-05-16  73  73.60548         47
#> 10 1952-11-14  73  73.10685         47
```

# Join [microorganisms](https://amr-for-r.org/release/3.0.1/reference/microorganisms.md) to a Data Set

Join the data set
[microorganisms](https://amr-for-r.org/release/3.0.1/reference/microorganisms.md)
easily to an existing data set or to a
[character](https://rdrr.io/r/base/character.html) vector.

## Usage

``` r
inner_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

left_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

right_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

full_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

semi_join_microorganisms(x, by = NULL, ...)

anti_join_microorganisms(x, by = NULL, ...)
```

## Arguments

- x:

  Existing data set to join, or
  [character](https://rdrr.io/r/base/character.html) vector. In case of
  a [character](https://rdrr.io/r/base/character.html) vector, the
  resulting [data.frame](https://rdrr.io/r/base/data.frame.html) will
  contain a column 'x' with these values.

- by:

  A variable to join by - if left empty will search for a column with
  class [`mo`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md)
  (created with
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md))
  or will be `"mo"` if that column name exists in `x`, could otherwise
  be a column name of `x` with values that exist in `microorganisms$mo`
  (such as `by = "bacteria_id"`), or another column in
  [microorganisms](https://amr-for-r.org/release/3.0.1/reference/microorganisms.md)
  (but then it should be named, like
  `by = c("bacteria_id" = "fullname")`).

- suffix:

  If there are non-joined duplicate variables in `x` and `y`, these
  suffixes will be added to the output to disambiguate them. Should be a
  [character](https://rdrr.io/r/base/character.html) vector of length 2.

- ...:

  Ignored, only in place to allow future extensions.

## Value

a [data.frame](https://rdrr.io/r/base/data.frame.html)

## Details

**Note:** As opposed to the `join()` functions of `dplyr`,
[character](https://rdrr.io/r/base/character.html) vectors are supported
and at default existing columns will get a suffix `"2"` and the newly
joined columns will not get a suffix.

If the `dplyr` package is installed, their join functions will be used.
Otherwise, the much slower
[`merge()`](https://rdrr.io/pkg/data.table/man/merge.html) and
[`interaction()`](https://rdrr.io/r/base/interaction.html) functions
from base R will be used.

## Examples

``` r
left_join_microorganisms(as.mo("K. pneumoniae"))
left_join_microorganisms("B_KLBSL_PNMN")

df <- data.frame(
  date = seq(
    from = as.Date("2018-01-01"),
    to = as.Date("2018-01-07"),
    by = 1
  ),
  bacteria = as.mo(c(
    "S. aureus", "MRSA", "MSSA", "STAAUR",
    "E. coli", "E. coli", "E. coli"
  )),
  stringsAsFactors = FALSE
)
colnames(df)

df_joined <- left_join_microorganisms(df, "bacteria")
colnames(df_joined)

# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    left_join_microorganisms() %>%
    colnames()
}
# }
```

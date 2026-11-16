# Check Availability of Columns

Easy check for data availability of all columns in a data set. This
makes it easy to get an idea of which antimicrobial combinations can be
used for calculation with e.g.
[`susceptibility()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md)
and
[`resistance()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md).

## Usage

``` r
availability(tbl, width = NULL)
```

## Arguments

- tbl:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) or
  [list](https://rdrr.io/r/base/list.html).

- width:

  Number of characters to present the visual availability - the default
  is filling the width of the console.

## Value

[data.frame](https://rdrr.io/r/base/data.frame.html) with column names
of `tbl` as row names

## Details

The function returns a
[data.frame](https://rdrr.io/r/base/data.frame.html) with columns
`"resistant"` and `"visual_resistance"`. The values in that columns are
calculated with
[`resistance()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md).

## Examples

``` r
availability(example_isolates)
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    filter(mo == as.mo("Escherichia coli")) %>%
    select_if(is.sir) %>%
    availability()
}
# }
```

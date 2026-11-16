# Filter Top *n* Microorganisms

This function filters a data set to include only the top *n*
microorganisms based on a specified property, such as taxonomic family
or genus. For example, it can filter a data set to the top 3 species, or
to any species in the top 5 genera, or to the top 3 species in each of
the top 5 genera.

## Usage

``` r
top_n_microorganisms(x, n, property = "species", n_for_each = NULL,
  col_mo = NULL, ...)
```

## Arguments

- x:

  A data frame containing microbial data.

- n:

  An integer specifying the maximum number of unique values of the
  `property` to include in the output.

- property:

  A character string indicating the microorganism property to use for
  filtering. Must be one of the column names of the
  [microorganisms](https://amr-for-r.org/release/3.0.1/reference/microorganisms.md)
  data set: "mo", "fullname", "status", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "subspecies", "rank", "ref",
  "oxygen_tolerance", "source", "lpsn", "lpsn_parent",
  "lpsn_renamed_to", "mycobank", "mycobank_parent",
  "mycobank_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to",
  "prevalence", or "snomed". If `NULL`, the raw values from `col_mo`
  will be used without transformation. When using `"species"` (default)
  or `"subpecies"`, the genus will be added to make sure each
  (sub)species still belongs to the right genus.

- n_for_each:

  An optional integer specifying the maximum number of rows to retain
  for each value of the selected property. If `NULL`, all rows within
  the top *n* groups will be included.

- col_mo:

  A character string indicating the column in `x` that contains
  microorganism names or codes. Defaults to the first column of class
  [`mo`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md).

- ...:

  Additional arguments passed on to
  [`mo_property()`](https://amr-for-r.org/release/3.0.1/reference/mo_property.md)
  when `property` is not `NULL`.

## Details

This function is useful for preprocessing data before creating
[antibiograms](https://amr-for-r.org/release/3.0.1/reference/antibiogram.md)
or other analyses that require focused subsets of microbial data. For
example, it can filter a data set to only include isolates from the top
10 species.

## See also

[`mo_property()`](https://amr-for-r.org/release/3.0.1/reference/mo_property.md),
[`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md),
[`antibiogram()`](https://amr-for-r.org/release/3.0.1/reference/antibiogram.md)

## Examples

``` r
# filter to the top 3 species:
top_n_microorganisms(example_isolates,
  n = 3
)

# filter to any species in the top 5 genera:
top_n_microorganisms(example_isolates,
  n = 5, property = "genus"
)

# filter to the top 3 species in each of the top 5 genera:
top_n_microorganisms(example_isolates,
  n = 5, property = "genus", n_for_each = 3
)
```

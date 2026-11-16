# Filter Top *n* Microorganisms

Filters a data set to include only the top *n* microorganisms based on a
specified property, such as taxonomic family or genus. For example, it
can filter a data set to the top 3 species, to any species in the top 5
genera, or to the top 3 species in each of the top 5 genera.

## Usage

``` r
top_n_microorganisms(
  x,
  n,
  property = "species",
  n_for_each = NULL,
  property_for_each = "species",
  col_mo = NULL,
  ...
)
```

## Arguments

- x:

  A data frame containing microbial data.

- n:

  A positive whole number specifying the maximum number of unique values
  of `property` to include in the output.

- property:

  A character string indicating the microorganism property to use for
  filtering. Must be one of the column names of the
  [microorganisms](https://amr-for-r.org/release/3.1.0/reference/microorganisms.md)
  data set: `"mo"`, `"fullname"`, `"status"`, `"domain"`, `"kingdom"`,
  `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"`, `"species"`,
  `"subspecies"`, `"rank"`, `"ref"`, `"oxygen_tolerance"`,
  `"morphology"`, `"source"`, `"lpsn"`, `"lpsn_parent"`,
  `"lpsn_renamed_to"`, `"mycobank"`, `"mycobank_parent"`,
  `"mycobank_renamed_to"`, `"gbif"`, `"gbif_parent"`,
  `"gbif_renamed_to"`, `"prevalence"`, or `"snomed"`. If `NULL`, the raw
  values from `col_mo` will be used without transformation. When using
  `"species"` (default) or `"subspecies"`, the genus is prepended to
  ensure each name is unambiguous.

- n_for_each:

  An optional positive whole number specifying the maximum number of
  distinct microorganism groups at the level of `property_for_each` to
  retain within each of the top *n* groups. Only used when
  `property_for_each` is also set.

- property_for_each:

  The microorganism property to use for sub-grouping within each top *n*
  group. Must be one of the column names of the
  [microorganisms](https://amr-for-r.org/release/3.1.0/reference/microorganisms.md)
  data set and at a strictly lower taxonomic rank than `property`
  (allowed order: domain \> kingdom \> phylum \> class \> order \>
  family \> genus \> species \> subspecies). Defaults to `"species"`.
  Only relevant when `n_for_each` is set.

- col_mo:

  A character string indicating the column in `x` that contains
  microorganism names or codes. Defaults to the first column of class
  [`mo`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md).

- ...:

  Additional arguments passed on to
  [`mo_property()`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md)
  when `property` is not `NULL`.

## Details

This function is useful for preprocessing data before creating
[antibiograms](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md)
or other analyses that require focused subsets of microbial data.

## See also

[`mo_property()`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md),
[`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md),
[`antibiogram()`](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md)

## Examples

``` r
# filter to the top 3 species:
top_n_microorganisms(example_isolates, n = 3)

# filter to any species in the top 5 genera:
top_n_microorganisms(example_isolates, n = 5, property = "genus")

# filter to the top 3 species in each of the top 5 genera:
top_n_microorganisms(example_isolates,
  n = 5, property = "genus", n_for_each = 3
)

# filter to the top 2 genera in each of the top 3 families:
top_n_microorganisms(example_isolates,
  n = 3, property = "family", n_for_each = 2, property_for_each = "genus"
)
```

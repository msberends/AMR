# Italicise Taxonomic Families, Genera, Species, Subspecies

According to the binomial nomenclature, the lowest four taxonomic levels
(family, genus, species, subspecies) should be printed in italics. This
function finds taxonomic names within strings and makes them italic.

## Usage

``` r
italicise_taxonomy(string, type = c("markdown", "ansi", "html"))

italicize_taxonomy(string, type = c("markdown", "ansi", "html"))
```

## Arguments

- string:

  A [character](https://rdrr.io/r/base/character.html) (vector).

- type:

  Type of conversion of the taxonomic names, either "markdown", "html"
  or "ansi", see *Details*.

## Details

This function finds the taxonomic names and makes them italic based on
the [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
data set.

The taxonomic names can be italicised using markdown (the default) by
adding `*` before and after the taxonomic names, or `<i>` and `</i>`
when using html. When using 'ansi', ANSI colours will be added using
`\033[3m` before and `\033[23m` after the taxonomic names. If multiple
ANSI colours are not available, no conversion will occur.

This function also supports abbreviation of the genus if it is followed
by a species, such as "E. coli" and "K. pneumoniae ozaenae".

## Examples

``` r
italicise_taxonomy("An overview of Staphylococcus aureus isolates")
#> [1] "An overview of *Staphylococcus aureus* isolates"
italicise_taxonomy("An overview of S. aureus isolates")
#> [1] "An overview of *S. aureus* isolates"

cat(italicise_taxonomy("An overview of S. aureus isolates", type = "ansi"))
#> An overview of S. aureus isolates
```

# Vectorised Pattern Matching with Keyboard Shortcut

Convenient wrapper around [`grepl()`](https://rdrr.io/r/base/grep.html)
to match a pattern: `x %like% pattern`. It always returns a
[`logical`](https://rdrr.io/r/base/logical.html) vector and is always
case-insensitive (use `x %like_case% pattern` for case-sensitive
matching). Also, `pattern` can be as long as `x` to compare items of
each index in both vectors, or they both can have the same length to
iterate over all cases.

## Usage

``` r
like(x, pattern, ignore.case = TRUE)

x %like% pattern

x %unlike% pattern

x %like_case% pattern

x %unlike_case% pattern
```

## Source

Idea from the [`like` function from the `data.table`
package](https://github.com/Rdatatable/data.table/blob/ec1259af1bf13fc0c96a1d3f9e84d55d8106a9a4/R/like.R),
although altered as explained in *Details*.

## Arguments

- x:

  A [character](https://rdrr.io/r/base/character.html) vector where
  matches are sought, or an object which can be coerced by
  [`as.character()`](https://rdrr.io/r/base/character.html) to a
  [character](https://rdrr.io/r/base/character.html) vector.

- pattern:

  A [character](https://rdrr.io/r/base/character.html) vector containing
  regular expressions (or a
  [character](https://rdrr.io/r/base/character.html) string for
  `fixed = TRUE`) to be matched in the given
  [character](https://rdrr.io/r/base/character.html) vector. Coerced by
  [`as.character()`](https://rdrr.io/r/base/character.html) to a
  [character](https://rdrr.io/r/base/character.html) string if possible.

- ignore.case:

  If `FALSE`, the pattern matching is *case sensitive* and if `TRUE`,
  case is ignored during matching.

## Value

A [logical](https://rdrr.io/r/base/logical.html) vector

## Details

These `like()` and `%like%`/`%unlike%` functions:

- Are case-insensitive (use `%like_case%`/`%unlike_case%` for
  case-sensitive matching)

- Support multiple patterns

- Check if `pattern` is a valid regular expression and sets
  `fixed = TRUE` if not, to greatly improve speed (vectorised over
  `pattern`)

- Always use compatibility with Perl unless `fixed = TRUE`, to greatly
  improve speed

Using RStudio? The `%like%`/`%unlike%` functions can also be directly
inserted in your code from the Addins menu and can have its own keyboard
shortcut like `Shift+Ctrl+L` or `Shift+Cmd+L` (see menu `Tools` \>
`Modify Keyboard Shortcuts...`). If you keep pressing your shortcut, the
inserted text will be iterated over `%like%` -\> `%unlike%` -\>
`%like_case%` -\> `%unlike_case%`.

## See also

[`grepl()`](https://rdrr.io/r/base/grep.html)

## Examples

``` r
# data.table has a more limited version of %like%, so unload it:
try(detach("package:data.table", unload = TRUE), silent = TRUE)
#> Warning: ‘data.table’ namespace cannot be unloaded:
#>   namespace ‘data.table’ is imported by ‘prodlim’ so cannot be unloaded

a <- "This is a test"
b <- "TEST"
a %like% b
#> [1] TRUE
b %like% a
#> [1] FALSE

# also supports multiple patterns
a <- c("Test case", "Something different", "Yet another thing")
b <- c("case", "diff", "yet")
a %like% b
#> [1] TRUE TRUE TRUE
a %unlike% b
#> [1] FALSE FALSE FALSE

a[1] %like% b
#> [1]  TRUE FALSE FALSE
a %like% b[1]
#> [1]  TRUE FALSE FALSE

# \donttest{
# get isolates whose name start with 'Entero' (case-insensitive)
example_isolates[which(mo_name() %like% "^entero"), ]
#> ℹ Using column 'mo' as input for `mo_name()`
#> # A tibble: 106 × 46
#>    date       patient   age gender ward    mo            PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>   <mo>          <sir> <sir> <sir> <sir>
#>  1 2002-02-21 4FC193     69 M      Clinic… B_ENTRC_FACM    NA    NA    NA    NA 
#>  2 2002-04-08 130252     78 M      ICU     B_ENTRC_FCLS    NA    NA    NA    NA 
#>  3 2002-06-23 798871     82 M      Clinic… B_ENTRC_FCLS    NA    NA    NA    NA 
#>  4 2002-06-23 798871     82 M      Clinic… B_ENTRC_FCLS    NA    NA    NA    NA 
#>  5 2003-04-20 6BC362     62 M      ICU     B_ENTRC         NA    NA    NA    NA 
#>  6 2003-04-21 6BC362     62 M      ICU     B_ENTRC         NA    NA    NA    NA 
#>  7 2003-08-13 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#>  8 2003-08-13 F35553     52 M      ICU     B_ENTRC_FCLS    NA    NA    NA    NA 
#>  9 2003-09-05 F35553     52 M      ICU     B_ENTRC         NA    NA    NA    NA 
#> 10 2003-09-05 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#> # ℹ 96 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

if (require("dplyr")) {
  example_isolates %>%
    filter(mo_name() %like% "^ent")
}
#> ℹ Using column 'mo' as input for `mo_name()`
#> # A tibble: 106 × 46
#>    date       patient   age gender ward    mo            PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>   <mo>          <sir> <sir> <sir> <sir>
#>  1 2002-02-21 4FC193     69 M      Clinic… B_ENTRC_FACM    NA    NA    NA    NA 
#>  2 2002-04-08 130252     78 M      ICU     B_ENTRC_FCLS    NA    NA    NA    NA 
#>  3 2002-06-23 798871     82 M      Clinic… B_ENTRC_FCLS    NA    NA    NA    NA 
#>  4 2002-06-23 798871     82 M      Clinic… B_ENTRC_FCLS    NA    NA    NA    NA 
#>  5 2003-04-20 6BC362     62 M      ICU     B_ENTRC         NA    NA    NA    NA 
#>  6 2003-04-21 6BC362     62 M      ICU     B_ENTRC         NA    NA    NA    NA 
#>  7 2003-08-13 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#>  8 2003-08-13 F35553     52 M      ICU     B_ENTRC_FCLS    NA    NA    NA    NA 
#>  9 2003-09-05 F35553     52 M      ICU     B_ENTRC         NA    NA    NA    NA 
#> 10 2003-09-05 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#> # ℹ 96 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
# }
```

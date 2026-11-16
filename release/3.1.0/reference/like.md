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
# data.table has a more limited version of %like%, so detach it
# (without unloading it, as data.table does not support being unloaded):
try(detach("package:data.table"), silent = TRUE)

a <- "This is a test"
b <- "TEST"
a %like% b
b %like% a

# also supports multiple patterns
a <- c("Test case", "Something different", "Yet another thing")
b <- c("case", "diff", "yet")
a %like% b
a %unlike% b

a[1] %like% b
a %like% b[1]

# \donttest{
# get isolates whose name start with 'Entero' (case-insensitive)
example_isolates[which(mo_name() %like% "^entero"), ]

if (require("dplyr")) {
  example_isolates %>%
    filter(mo_name() %like% "^ent")
}
# }
```

# Retrieve Antimicrobial Drug Names and Doses from Clinical Text

Use this function on e.g. clinical texts from health care records. It
returns a [list](https://rdrr.io/r/base/list.html) with all
antimicrobial drugs, doses and forms of administration found in the
texts.

## Usage

``` r
ab_from_text(text, type = c("drug", "dose", "administration"),
  collapse = NULL, translate_ab = FALSE, thorough_search = NULL,
  info = interactive(), ...)
```

## Arguments

- text:

  Text to analyse.

- type:

  Type of property to search for, either `"drug"`, `"dose"` or
  `"administration"`, see *Examples*.

- collapse:

  A [character](https://rdrr.io/r/base/character.html) to pass on to
  `paste(, collapse = ...)` to only return one
  [character](https://rdrr.io/r/base/character.html) per element of
  `text`, see *Examples*.

- translate_ab:

  If `type = "drug"`: a column name of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set to translate the antibiotic abbreviations to, using
  [`ab_property()`](https://amr-for-r.org/reference/ab_property.md). The
  default is `FALSE`. Using `TRUE` is equal to using "name".

- thorough_search:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the input must be extensively searched for misspelling and other
  faulty input values. Setting this to `TRUE` will take considerably
  more time than when using `FALSE`. At default, it will turn `TRUE`
  when all input elements contain a maximum of three words.

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  progress bar should be printed - the default is `TRUE` only in
  interactive mode.

- ...:

  Arguments passed on to
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

## Value

A [list](https://rdrr.io/r/base/list.html), or a
[character](https://rdrr.io/r/base/character.html) if `collapse` is not
`NULL`

## Details

This function is also internally used by
[`as.ab()`](https://amr-for-r.org/reference/as.ab.md), although it then
only searches for the first drug name and will throw a note if more drug
names could have been returned. Note: the
[`as.ab()`](https://amr-for-r.org/reference/as.ab.md) function may use
very long regular expression to match brand names of antimicrobial
drugs. This may fail on some systems.

### Argument `type`

At default, the function will search for antimicrobial drug names. All
text elements will be searched for official names, ATC codes and brand
names. As it uses [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
internally, it will correct for misspelling.

With `type = "dose"` (or similar, like "dosing", "doses"), all text
elements will be searched for
[numeric](https://rdrr.io/r/base/numeric.html) values that are higher
than 100 and do not resemble years. The output will be
[numeric](https://rdrr.io/r/base/numeric.html). It supports any unit (g,
mg, IE, etc.) and multiple values in one clinical text, see *Examples*.

With `type = "administration"` (or abbreviations, like "admin", "adm"),
all text elements will be searched for a form of drug administration. It
supports the following forms (including common abbreviations): buccal,
implant, inhalation, instillation, intravenous, nasal, oral, parenteral,
rectal, sublingual, transdermal and vaginal. Abbreviations for oral
(such as 'po', 'per os') will become "oral", all values for intravenous
(such as 'iv', 'intraven') will become "iv". It supports multiple values
in one clinical text, see *Examples*.

### Argument `collapse`

Without using `collapse`, this function will return a
[list](https://rdrr.io/r/base/list.html). This can be convenient to use
e.g. inside a
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)):  
`df %>% mutate(abx = ab_from_text(clinical_text))`

The returned AB codes can be transformed to official names, groups, etc.
with all [`ab_*`](https://amr-for-r.org/reference/ab_property.md)
functions such as
[`ab_name()`](https://amr-for-r.org/reference/ab_property.md) and
[`ab_group()`](https://amr-for-r.org/reference/ab_property.md), or by
using the `translate_ab` argument.

With using `collapse`, this function will return a
[character](https://rdrr.io/r/base/character.html):  
`df %>% mutate(abx = ab_from_text(clinical_text, collapse = "|"))`

## Examples

``` r
# mind the bad spelling of amoxicillin in this line,
# straight from a true health care record:
ab_from_text("28/03/2020 regular amoxicilliin 500mg po tid")
#> [[1]]
#> Class 'ab'
#> [1] AMX
#> 

ab_from_text("500 mg amoxi po and 400mg cipro iv")
#> [[1]]
#> Class 'ab'
#> [1] AMX CIP
#> 
ab_from_text("500 mg amoxi po and 400mg cipro iv", type = "dose")
#> [[1]]
#> [1] 500 400
#> 
ab_from_text("500 mg amoxi po and 400mg cipro iv", type = "admin")
#> [[1]]
#> [1] "oral" "iv"  
#> 

ab_from_text("500 mg amoxi po and 400mg cipro iv", collapse = ", ")
#> [1] "AMX, CIP"
# \donttest{
# if you want to know which antibiotic groups were administered, do e.g.:
abx <- ab_from_text("500 mg amoxi po and 400mg cipro iv")
ab_group(abx[[1]])
#> [1] "Aminopenicillins" "Fluoroquinolones"

if (require("dplyr")) {
  tibble(clinical_text = c(
    "given 400mg cipro and 500 mg amox",
    "started on doxy iv today"
  )) %>%
    mutate(
      abx_codes = ab_from_text(clinical_text),
      abx_doses = ab_from_text(clinical_text, type = "doses"),
      abx_admin = ab_from_text(clinical_text, type = "admin"),
      abx_coll = ab_from_text(clinical_text, collapse = "|"),
      abx_coll_names = ab_from_text(clinical_text,
        collapse = "|",
        translate_ab = "name"
      ),
      abx_coll_doses = ab_from_text(clinical_text,
        type = "doses",
        collapse = "|"
      ),
      abx_coll_admin = ab_from_text(clinical_text,
        type = "admin",
        collapse = "|"
      )
    )
}
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
#> # A tibble: 2 × 8
#>   clinical_text            abx_codes abx_doses abx_admin abx_coll abx_coll_names
#>   <chr>                    <list>    <list>    <list>    <chr>    <chr>         
#> 1 given 400mg cipro and 5… <ab [2]>  <dbl [2]> <chr [1]> CIP|AMX  Ciprofloxacin…
#> 2 started on doxy iv today <ab [1]>  <dbl [1]> <chr [1]> DOX      Doxycycline   
#> # ℹ 2 more variables: abx_coll_doses <chr>, abx_coll_admin <chr>
# }
```

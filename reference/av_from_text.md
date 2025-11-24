# Retrieve Antiviral Drug Names and Doses from Clinical Text

Use this function on e.g. clinical texts from health care records. It
returns a [list](https://rdrr.io/r/base/list.html) with all antiviral
drugs, doses and forms of administration found in the texts.

## Usage

``` r
av_from_text(text, type = c("drug", "dose", "administration"),
  collapse = NULL, translate_av = FALSE, thorough_search = NULL,
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

- translate_av:

  If `type = "drug"`: a column name of the
  [antivirals](https://amr-for-r.org/reference/antimicrobials.md) data
  set to translate the antibiotic abbreviations to, using
  [`av_property()`](https://amr-for-r.org/reference/av_property.md). The
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
  [`as.av()`](https://amr-for-r.org/reference/as.av.md).

## Value

A [list](https://rdrr.io/r/base/list.html), or a
[character](https://rdrr.io/r/base/character.html) if `collapse` is not
`NULL`

## Details

This function is also internally used by
[`as.av()`](https://amr-for-r.org/reference/as.av.md), although it then
only searches for the first drug name and will throw a note if more drug
names could have been returned. Note: the
[`as.av()`](https://amr-for-r.org/reference/as.av.md) function may use
very long regular expression to match brand names of antiviral drugs.
This may fail on some systems.

### Argument `type`

At default, the function will search for antiviral drug names. All text
elements will be searched for official names, ATC codes and brand names.
As it uses [`as.av()`](https://amr-for-r.org/reference/as.av.md)
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
`df %>% mutate(avx = av_from_text(clinical_text))`

The returned AV codes can be transformed to official names, groups, etc.
with all [`av_*`](https://amr-for-r.org/reference/av_property.md)
functions such as
[`av_name()`](https://amr-for-r.org/reference/av_property.md) and
[`av_group()`](https://amr-for-r.org/reference/av_property.md), or by
using the `translate_av` argument.

With using `collapse`, this function will return a
[character](https://rdrr.io/r/base/character.html):  
`df %>% mutate(avx = av_from_text(clinical_text, collapse = "|"))`

## Examples

``` r
av_from_text("28/03/2020 valaciclovir po tid")
#> [[1]]
#> Class 'av'
#> [1] VALA
#> 
av_from_text("28/03/2020 valaciclovir po tid", type = "admin")
#> [[1]]
#> [1] "oral"
#> 
```

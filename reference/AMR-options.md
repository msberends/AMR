# Options for the AMR package

This is an overview of all the package-specific
[`options()`](https://rdrr.io/r/base/options.html) you can set in the
`AMR` package.

## Options

- `AMR_antibiogram_formatting_type`  
  A [numeric](https://rdrr.io/r/base/numeric.html) (1-22) to use in
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md), to
  indicate which formatting type to use.

- `AMR_breakpoint_type`  
  A [character](https://rdrr.io/r/base/character.html) to use in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md), to indicate
  which breakpoint type to use. This must be either "ECOFF", "animal",
  or "human".

- `AMR_capped_mic_handling`  
  A [character](https://rdrr.io/r/base/character.html) to use in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md), to indicate
  how capped MIC values (`<`, `<=`, `>`, `>=`) should be interpreted.
  Must be one of `"none"`, `"conservative"`, `"standard"`, or
  `"lenient"` - the default is `"conservative"`.

- `AMR_cleaning_regex`  
  A [regular expression](https://rdrr.io/r/base/regex.html)
  (case-insensitive) to use in
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all
  [`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions, to
  clean the user input. The default is the outcome of
  [`mo_cleaning_regex()`](https://amr-for-r.org/reference/as.mo.md),
  which removes texts between brackets and texts such as "species" and
  "serovar".

- `AMR_custom_ab`  
  A file location to an RDS file, to use custom antimicrobial drugs with
  this package. This is explained in
  [`add_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md).

- `AMR_custom_mo`  
  A file location to an RDS file, to use custom microorganisms with this
  package. This is explained in
  [`add_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md).

- `AMR_eucastrules`  
  A [character](https://rdrr.io/r/base/character.html) to set the
  default types of rules for
  [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
  function, must be one or more of: `"breakpoints"`, `"expert"`,
  `"other"`, `"custom"`, `"all"`, and defaults to
  `c("breakpoints", "expert")`.

- `AMR_guideline`  
  A [character](https://rdrr.io/r/base/character.html) to set the
  default guideline for interpreting MIC values and disk diffusion
  diameters with
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md). Can be only
  the guideline name (e.g., `"CLSI"`) or the name with a year (e.g.
  `"CLSI 2019"`). The default to the latest implemented EUCAST
  guideline, currently `"EUCAST 2025"`. Supported guideline are
  currently EUCAST (2011-2025) and CLSI (2011-2025).

- `AMR_ignore_pattern`  
  A [regular expression](https://rdrr.io/r/base/regex.html) to ignore
  (i.e., make `NA`) any match given in
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all
  [`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions.

- `AMR_include_PKPD`  
  A [logical](https://rdrr.io/r/base/logical.html) to use in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md), to indicate
  that PK/PD clinical breakpoints must be applied as a last resort - the
  default is `TRUE`.

- `AMR_substitute_missing_r_breakpoint`  
  A [logical](https://rdrr.io/r/base/logical.html) to use in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md), to indicate
  that missing R breakpoints must be substituted with `"R"` - the
  default is `FALSE`.

- `AMR_include_screening`  
  A [logical](https://rdrr.io/r/base/logical.html) to use in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md), to indicate
  that clinical breakpoints for screening are allowed - the default is
  `FALSE`.

- `AMR_keep_synonyms`  
  A [logical](https://rdrr.io/r/base/logical.html) to use in
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all
  [`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions, to
  indicate if old, previously valid taxonomic names must be preserved
  and not be corrected to currently accepted names. The default is
  `FALSE`.

- `AMR_locale`  
  A [character](https://rdrr.io/r/base/character.html) to set the
  language for the `AMR` package, can be one of these supported language
  names or [ISO 639-1 codes](https://en.wikipedia.org/wiki/ISO_639-1):
  English (en), Arabic (ar), Bengali (bn), Chinese (zh), Czech (cs),
  Danish (da), Dutch (nl), Finnish (fi), French (fr), German (de), Greek
  (el), Hindi (hi), Indonesian (id), Italian (it), Japanese (ja), Korean
  (ko), Norwegian (no), Polish (pl), Portuguese (pt), Romanian (ro),
  Russian (ru), Spanish (es), Swahili (sw), Swedish (sv), Turkish (tr),
  Ukrainian (uk), Urdu (ur), or Vietnamese (vi). The default is the
  current system language (if supported, English otherwise).

- `AMR_mo_source`  
  A file location for a manual code list to be used in
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all
  [`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions.
  This is explained in
  [`set_mo_source()`](https://amr-for-r.org/reference/mo_source.md).

## Saving Settings Between Sessions

Settings in R are not saved globally and are thus lost when R is exited.
You can save your options to your own `.Rprofile` file, which is a
user-specific file. You can edit it using:

      utils::file.edit("~/.Rprofile")

In this file, you can set options such as...

     options(AMR_locale = "pt")
     options(AMR_include_PKPD = TRUE)

...to add Portuguese language support of antimicrobials, and allow PK/PD
rules when interpreting MIC values with
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md).

### Share Options Within Team

For a more global approach, e.g. within a (data) team, save an options
file to a remote file location, such as a shared network drive, and have
each user read in this file automatically at start-up. This would work
in this way:

1.  Save a plain text file to e.g. "X:/team_folder/R_options.R" and fill
    it with preferred settings.

2.  For each user, open the `.Rprofile` file using
    `utils::file.edit("~/.Rprofile")` and put in there:

          source("X:/team_folder/R_options.R")

3.  Reload R/RStudio and check the settings with
    [`getOption()`](https://rdrr.io/r/base/options.html), e.g.
    `getOption("AMR_locale")` if you have set that value.

Now the team settings are configured in only one place, and can be
maintained there.

# Translate Strings from the AMR Package

For language-dependent output of `AMR` functions, such as
[`mo_name()`](https://amr-for-r.org/reference/mo_property.md),
[`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md),
[`mo_type()`](https://amr-for-r.org/reference/mo_property.md) and
[`ab_name()`](https://amr-for-r.org/reference/ab_property.md).

## Usage

``` r
get_AMR_locale()

set_AMR_locale(language)

reset_AMR_locale()

translate_AMR(x, language = get_AMR_locale())
```

## Arguments

- language:

  Language to choose. Use one of these supported language names or [ISO
  639-1 codes](https://en.wikipedia.org/wiki/ISO_639-1): English (en),
  Arabic (ar), Bengali (bn), Chinese (zh), Czech (cs), Danish (da),
  Dutch (nl), Finnish (fi), French (fr), German (de), Greek (el), Hindi
  (hi), Indonesian (id), Italian (it), Japanese (ja), Korean (ko),
  Norwegian (no), Polish (pl), Portuguese (pt), Romanian (ro), Russian
  (ru), Spanish (es), Swahili (sw), Swedish (sv), Turkish (tr),
  Ukrainian (uk), Urdu (ur), or Vietnamese (vi).

- x:

  Text to translate.

## Details

The currently 28 supported languages are English (en), Arabic (ar),
Bengali (bn), Chinese (zh), Czech (cs), Danish (da), Dutch (nl), Finnish
(fi), French (fr), German (de), Greek (el), Hindi (hi), Indonesian (id),
Italian (it), Japanese (ja), Korean (ko), Norwegian (no), Polish (pl),
Portuguese (pt), Romanian (ro), Russian (ru), Spanish (es), Swahili
(sw), Swedish (sv), Turkish (tr), Ukrainian (uk), Urdu (ur), and
Vietnamese (vi). All these languages have translations available for all
antimicrobial drugs and colloquial microorganism names.

To permanently silence the once-per-session language note on a
non-English operating system, you can set the package option
[`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md) in your
`.Rprofile` file like this:

    # Open .Rprofile file
    utils::file.edit("~/.Rprofile")

    # Then add e.g. Italian support to that file using:
    options(AMR_locale = "Italian")

And then save the file.

Please read about adding or updating a language in [our
Wiki](https://github.com/msberends/AMR/wiki/).

### Changing the Default Language

The system language will be used at default (as returned by
`Sys.getenv("LANG")` or, if `LANG` is not set,
[`Sys.getlocale("LC_COLLATE")`](https://rdrr.io/r/base/locales.html)),
if that language is supported. But the language to be used can be
overwritten in two ways and will be checked in this order:

1.  Setting the package option
    [`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md),
    either by using e.g. `set_AMR_locale("German")` or by running e.g.
    `options(AMR_locale = "German")`.

    Note that setting an R option only works in the same session. Save
    the command `options(AMR_locale = "(your language)")` to your
    `.Rprofile` file to apply it for every session. Run
    `utils::file.edit("~/.Rprofile")` to edit your `.Rprofile` file.

2.  Setting the system variable `LANGUAGE` or `LANG`, e.g. by adding
    `LANGUAGE="de_DE.utf8"` to your `.Renviron` file in your home
    directory.

Thus, if the package option
[`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md) is set,
the system variables `LANGUAGE` and `LANG` will be ignored.

## Examples

``` r
# Current settings (based on system language)
ab_name("Ciprofloxacin")
#> [1] "Ciprofloxacin"
mo_name("Coagulase-negative Staphylococcus (CoNS)")
#> [1] "Coagulase-negative Staphylococcus (CoNS)"

# setting another language
set_AMR_locale("Dutch")
#> ℹ Using Dutch (Nederlands) for the AMR package for this session.
ab_name("Ciprofloxacin")
#> [1] "Ciprofloxacine"
mo_name("Coagulase-negative Staphylococcus (CoNS)")
#> [1] "Coagulase-negatieve Staphylococcus (CNS)"

# setting yet another language
set_AMR_locale("German")
#> ℹ Using German (Deutsch) for the AMR package for this session.
ab_name("Ciprofloxacin")
#> [1] "Ciprofloxacin"
mo_name("Coagulase-negative Staphylococcus (CoNS)")
#> [1] "Koagulase-negative Staphylococcus (KNS)"

# set_AMR_locale() understands endonyms, English exonyms, and ISO 639-1:
set_AMR_locale("Deutsch")
#> ℹ Using German (Deutsch) for the AMR package for this session.
set_AMR_locale("German")
#> ℹ Using German (Deutsch) for the AMR package for this session.
set_AMR_locale("de")
#> ℹ Using German (Deutsch) for the AMR package for this session.
ab_name("amox/clav")
#> [1] "Amoxicillin/Clavulansäure"

# reset to system default
reset_AMR_locale()
#> ℹ Using the English language (English) for the AMR package for this
#>   session.
ab_name("amox/clav")
#> [1] "Amoxicillin/clavulanic acid"
```

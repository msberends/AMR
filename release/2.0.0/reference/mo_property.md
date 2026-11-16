# Get Properties of a Microorganism

Use these functions to return a specific property of a microorganism
based on the latest accepted taxonomy. All input values will be
evaluated internally with
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md),
which makes it possible to use microbial abbreviations, codes and names
as input. See *Examples*.

## Usage

``` r
mo_name(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_fullname(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_shortname(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_subspecies(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_species(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_genus(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_family(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_order(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_class(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_phylum(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_kingdom(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_domain(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_type(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_status(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_pathogenicity(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_gramstain(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_is_gram_negative(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_is_gram_positive(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_is_yeast(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_is_intrinsic_resistant(
  x,
  ab,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_snomed(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_ref(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_authors(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_year(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_lpsn(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_gbif(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_rank(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_taxonomy(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_synonyms(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_current(x, language = get_AMR_locale(), ...)

mo_info(
  x,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_url(
  x,
  open = FALSE,
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)

mo_property(
  x,
  property = "fullname",
  language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  ...
)
```

## Arguments

- x:

  any [character](https://rdrr.io/r/base/character.html) (vector) that
  can be coerced to a valid microorganism code with
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md).
  Can be left blank for auto-guessing the column containing
  microorganism codes if used in a data set, see *Examples*.

- language:

  language to translate text like "no growth", which defaults to the
  system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/2.0.0/reference/translate.md))

- keep_synonyms:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate if old,
  previously valid taxonomic names must be preserved and not be
  corrected to currently accepted names. The default is `FALSE`, which
  will return a note if old taxonomic names were processed. The default
  can be set with the [package
  option](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md)
  [`AMR_keep_synonyms`](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md),
  i.e. `options(AMR_keep_synonyms = TRUE)` or
  `options(AMR_keep_synonyms = FALSE)`.

- ...:

  other arguments passed on to
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md),
  such as 'minimum_matching_score', 'ignore_pattern', and
  'remove_from_input'

- ab:

  any (vector of) text that can be coerced to a valid antibiotic drug
  code with
  [`as.ab()`](https://amr-for-r.org/release/2.0.0/reference/as.ab.md)

- open:

  browse the URL using
  [`browseURL()`](https://rdrr.io/r/utils/browseURL.html)

- property:

  one of the column names of the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set: "mo", "fullname", "status", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "subspecies", "rank", "ref",
  "source", "lpsn", "lpsn_parent", "lpsn_renamed_to", "gbif",
  "gbif_parent", "gbif_renamed_to", "prevalence", or "snomed", or must
  be `"shortname"`

## Value

- An [integer](https://rdrr.io/r/base/integer.html) in case of
  `mo_year()`

- An [ordered factor](https://rdrr.io/pkg/data.table/man/fctr.html) in
  case of `mo_pathogenicity()`

- A [list](https://rdrr.io/r/base/list.html) in case of `mo_taxonomy()`,
  `mo_synonyms()`, `mo_snomed()` and `mo_info()`

- A named [character](https://rdrr.io/r/base/character.html) in case of
  `mo_url()`

- A [character](https://rdrr.io/r/base/character.html) in all other
  cases

## Details

All functions will, at default, **not** keep old taxonomic properties,
as synonyms are automatically replaced with the current taxonomy. Take
for example *Enterobacter aerogenes*, which was initially named in 1960
but renamed to *Klebsiella aerogenes* in 2017:

- `mo_genus("Enterobacter aerogenes")` will return `"Klebsiella"` (with
  a note about the renaming)

- `mo_genus("Enterobacter aerogenes", keep_synonyms = TRUE)` will return
  `"Enterobacter"` (with a once-per-session warning that the name is
  outdated)

- `mo_ref("Enterobacter aerogenes")` will return
  `"Tindall et al., 2017"` (with a note)

- `mo_ref("Enterobacter aerogenes", keep_synonyms = TRUE)` will return
  `"Hormaeche et al., 1960"` (with a warning)

The short name (`mo_shortname()`) returns the first character of the
genus and the full species, such as `"E. coli"`, for species and
subspecies. Exceptions are abbreviations of staphylococci (such as
*"CoNS"*, Coagulase-Negative Staphylococci) and beta-haemolytic
streptococci (such as *"GBS"*, Group B Streptococci). Please bear in
mind that e.g. *E. coli* could mean *Escherichia coli* (kingdom of
Bacteria) as well as *Entamoeba coli* (kingdom of Protozoa). Returning
to the full name will be done using
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)
internally, giving priority to bacteria and human pathogens, i.e.
`"E. coli"` will be considered *Escherichia coli*. As a result,
`mo_fullname(mo_shortname("Entamoeba coli"))` returns
`"Escherichia coli"`.

Since the top-level of the taxonomy is sometimes referred to as
'kingdom' and sometimes as 'domain', the functions `mo_kingdom()` and
`mo_domain()` return the exact same results.

Determination of human pathogenicity (`mo_pathogenicity()`) is strongly
based on Bartlett *et al.* (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ). This
function returns a
[factor](https://rdrr.io/pkg/data.table/man/fctr.html) with the levels
*Pathogenic*, *Potentially pathogenic*, *Non-pathogenic*, and *Unknown*.

Determination of the Gram stain (`mo_gramstain()`) will be based on the
taxonomic kingdom and phylum. Originally, Cavalier-Smith defined the
so-called subkingdoms Negibacteria and Posibacteria (2002, [PMID
11837318](https://pubmed.ncbi.nlm.nih.gov/11837318/)), and only
considered these phyla as Posibacteria: Actinobacteria, Chloroflexi,
Firmicutes, and Tenericutes. These phyla were later renamed to
Actinomycetota, Chloroflexota, Bacillota, and Mycoplasmatota (2021,
[PMID 34694987](https://pubmed.ncbi.nlm.nih.gov/34694987/)). Bacteria in
these phyla are considered Gram-positive in this `AMR` package, except
for members of the class Negativicutes (within phylum Bacillota) which
are Gram-negative. All other bacteria are considered Gram-negative.
Species outside the kingdom of Bacteria will return a value `NA`.
Functions `mo_is_gram_negative()` and `mo_is_gram_positive()` always
return `TRUE` or `FALSE` (or `NA` when the input is `NA` or the MO code
is `UNKNOWN`), thus always return `FALSE` for species outside the
taxonomic kingdom of Bacteria.

Determination of yeasts (`mo_is_yeast()`) will be based on the taxonomic
kingdom and class. *Budding yeasts* are fungi of the phylum Ascomycota,
class Saccharomycetes (also called Hemiascomycetes). *True yeasts* are
aggregated into the underlying order Saccharomycetales. Thus, for all
microorganisms that are member of the taxonomic class Saccharomycetes,
the function will return `TRUE`. It returns `FALSE` otherwise (or `NA`
when the input is `NA` or the MO code is `UNKNOWN`).

Determination of intrinsic resistance (`mo_is_intrinsic_resistant()`)
will be based on the
[intrinsic_resistant](https://amr-for-r.org/release/2.0.0/reference/intrinsic_resistant.md)
data set, which is based on ['EUCAST Expert Rules' and 'EUCAST Intrinsic
Resistance and Unusual Phenotypes'
v3.3](https://www.eucast.org/expert_rules_and_expected_phenotypes/)
(2021). The `mo_is_intrinsic_resistant()` function can be vectorised
over both argument `x` (input for microorganisms) and `ab` (input for
antibiotics).

The function `mo_url()` will return the direct URL to the online
database entry, which also shows the scientific reference of the
concerned species.

SNOMED codes (`mo_snomed()`) are from the version of 1 July, 2021. See
*Source* and the
[microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
data set for more info.

Old taxonomic names (so-called 'synonyms') can be retrieved with
`mo_synonyms()` (which will have the scientific reference as
[name](https://rdrr.io/r/base/names.html)), the current taxonomic name
can be retrieved with `mo_current()`. Both functions return full names.

All output [will be
translated](https://amr-for-r.org/release/2.0.0/reference/translate.md)
where possible.

## Matching Score for Microorganisms

This function uses
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)
internally, which uses an advanced algorithm to translate arbitrary user
input to valid taxonomy using a so-called matching score. You can read
about this public algorithm on the [MO matching score
page](https://amr-for-r.org/release/2.0.0/reference/mo_matching_score.md).

## Source

1.  Berends MS *et al.* (2022). **AMR: An R Package for Working with
    Antimicrobial Resistance Data**. *Journal of Statistical Software*,
    104(3), 1-31;
    [doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)

2.  Becker K *et al.* (2014). **Coagulase-Negative Staphylococci.**
    *Clin Microbiol Rev.* 27(4): 870-926;
    [doi:10.1128/CMR.00109-13](https://doi.org/10.1128/CMR.00109-13)

3.  Becker K *et al.* (2019). **Implications of identifying the recently
    defined members of the *S. aureus* complex, *S. argenteus* and *S.
    schweitzeri*: A position paper of members of the ESCMID Study Group
    for staphylococci and Staphylococcal Diseases (ESGS).** *Clin
    Microbiol Infect*;
    [doi:10.1016/j.cmi.2019.02.028](https://doi.org/10.1016/j.cmi.2019.02.028)

4.  Becker K *et al.* (2020). **Emergence of coagulase-negative
    staphylococci** *Expert Rev Anti Infect Ther.* 18(4):349-366;
    [doi:10.1080/14787210.2020.1730813](https://doi.org/10.1080/14787210.2020.1730813)

5.  Lancefield RC (1933). **A serological differentiation of human and
    other groups of hemolytic streptococci**. *J Exp Med.* 57(4):
    571-95;
    [doi:10.1084/jem.57.4.571](https://doi.org/10.1084/jem.57.4.571)

6.  Berends MS *et al.* (2022). **Trends in Occurrence and Phenotypic
    Resistance of Coagulase-Negative Staphylococci (CoNS) Found in Human
    Blood in the Northern Netherlands between 2013 and 2019**
    *Microorganisms* 10(9), 1801;
    [doi:10.3390/microorganisms10091801](https://doi.org/10.3390/microorganisms10091801)

7.  Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing
    in Nomenclature (LPSN) moves to the DSMZ.** International Journal of
    Systematic and Evolutionary Microbiology, 70, 5607-5612;
    [doi:10.1099/ijsem.0.004332](https://doi.org/10.1099/ijsem.0.004332)
    . Accessed from <https://lpsn.dsmz.de> on 11 December, 2022.

8.  GBIF Secretariat (2022). GBIF Backbone Taxonomy. Checklist dataset
    [doi:10.15468/39omei](https://doi.org/10.15468/39omei) . Accessed
    from <https://www.gbif.org> on 11 December, 2022.

9.  Public Health Information Network Vocabulary Access and Distribution
    System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
    Value Set Name 'Microoganism', OID 2.16.840.1.114222.4.11.1009
    (v12). URL: <https://phinvads.cdc.gov>

10. Bartlett A *et al.* (2022). **A comprehensive list of bacterial
    pathogens infecting humans** *Microbiology* 168:001269;
    [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## Reference Data Publicly Available

All data sets in this `AMR` package (about microorganisms, antibiotics,
SIR interpretation, EUCAST rules, etc.) are publicly and freely
available for download in the following formats: R, MS Excel, Apache
Feather, Apache Parquet, SPSS, SAS, and Stata. We also provide
tab-separated plain text files that are machine-readable and suitable
for input in any software program, such as laboratory information
systems. Please visit [our website for the download
links](https://msberends.github.io/AMR/articles/datasets.html). The
actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## See also

Data set
[microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)

## Examples

``` r
# taxonomic tree -----------------------------------------------------------

mo_kingdom("Klebsiella pneumoniae")
mo_phylum("Klebsiella pneumoniae")
mo_class("Klebsiella pneumoniae")
mo_order("Klebsiella pneumoniae")
mo_family("Klebsiella pneumoniae")
mo_genus("Klebsiella pneumoniae")
mo_species("Klebsiella pneumoniae")
mo_subspecies("Klebsiella pneumoniae")


# full names and short names -----------------------------------------------

mo_name("Klebsiella pneumoniae")
mo_fullname("Klebsiella pneumoniae")
mo_shortname("Klebsiella pneumoniae")


# other properties ---------------------------------------------------------

mo_pathogenicity("Klebsiella pneumoniae")
mo_gramstain("Klebsiella pneumoniae")
mo_snomed("Klebsiella pneumoniae")
mo_type("Klebsiella pneumoniae")
mo_rank("Klebsiella pneumoniae")
mo_url("Klebsiella pneumoniae")
mo_is_yeast(c("Candida", "Trichophyton", "Klebsiella"))


# scientific reference -----------------------------------------------------

mo_ref("Klebsiella aerogenes")
mo_authors("Klebsiella aerogenes")
mo_year("Klebsiella aerogenes")
mo_lpsn("Klebsiella aerogenes")
mo_gbif("Klebsiella aerogenes")
mo_synonyms("Klebsiella aerogenes")


# abbreviations known in the field -----------------------------------------

mo_genus("MRSA")
mo_species("MRSA")
mo_shortname("VISA")
mo_gramstain("VISA")

mo_genus("EHEC")
mo_species("EIEC")
mo_name("UPEC")


# known subspecies ---------------------------------------------------------

mo_fullname("K. pneu rh")
mo_shortname("K. pneu rh")

# \donttest{
# Becker classification, see ?as.mo ----------------------------------------

mo_fullname("Staph epidermidis")
mo_fullname("Staph epidermidis", Becker = TRUE)
mo_shortname("Staph epidermidis")
mo_shortname("Staph epidermidis", Becker = TRUE)


# Lancefield classification, see ?as.mo ------------------------------------

mo_fullname("Strep agalactiae")
mo_fullname("Strep agalactiae", Lancefield = TRUE)
mo_shortname("Strep agalactiae")
mo_shortname("Strep agalactiae", Lancefield = TRUE)


# language support  --------------------------------------------------------

mo_gramstain("Klebsiella pneumoniae", language = "de") # German
mo_gramstain("Klebsiella pneumoniae", language = "nl") # Dutch
mo_gramstain("Klebsiella pneumoniae", language = "es") # Spanish
mo_gramstain("Klebsiella pneumoniae", language = "el") # Greek
mo_gramstain("Klebsiella pneumoniae", language = "uk") # Ukrainian

# mo_type is equal to mo_kingdom, but mo_kingdom will remain untranslated
mo_kingdom("Klebsiella pneumoniae")
mo_type("Klebsiella pneumoniae")
mo_kingdom("Klebsiella pneumoniae", language = "zh") # Chinese, no effect
mo_type("Klebsiella pneumoniae", language = "zh") # Chinese, translated

mo_fullname("S. pyogenes", Lancefield = TRUE, language = "de")
mo_fullname("S. pyogenes", Lancefield = TRUE, language = "uk")


# other --------------------------------------------------------------------

# gram stains and intrinsic resistance can be used as a filter in dplyr verbs
if (require("dplyr")) {
  example_isolates %>%
    filter(mo_is_gram_positive()) %>%
    count(mo_genus(), sort = TRUE)
}
if (require("dplyr")) {
  example_isolates %>%
    filter(mo_is_intrinsic_resistant(ab = "vanco")) %>%
    count(mo_genus(), sort = TRUE)
}

# get a list with the complete taxonomy (from kingdom to subspecies)
mo_taxonomy("Klebsiella pneumoniae")

# get a list with the taxonomy, the authors, Gram-stain,
# SNOMED codes, and URL to the online database
mo_info("Klebsiella pneumoniae")
# }
```

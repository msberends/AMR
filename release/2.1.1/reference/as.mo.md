# Transform Arbitrary Input to Valid Microbial Taxonomy

Use this function to get a valid microorganism code (`mo`) based on
arbitrary user input. Determination is done using intelligent rules and
the complete taxonomic tree of the kingdoms Animalia, Archaea, Bacteria,
and Protozoa, and most microbial species from the kingdom Fungi (see
*Source*). The input can be almost anything: a full name (like
`"Staphylococcus aureus"`), an abbreviated name (such as `"S. aureus"`),
an abbreviation known in the field (such as `"MRSA"`), or just a genus.
See *Examples*.

## Usage

``` r
as.mo(
  x,
  Becker = FALSE,
  Lancefield = FALSE,
  minimum_matching_score = NULL,
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  reference_df = get_mo_source(),
  ignore_pattern = getOption("AMR_ignore_pattern", NULL),
  cleaning_regex = getOption("AMR_cleaning_regex", mo_cleaning_regex()),
  language = get_AMR_locale(),
  info = interactive(),
  ...
)

is.mo(x)

mo_uncertainties()

mo_renamed()

mo_failures()

mo_reset_session()

mo_cleaning_regex()
```

## Arguments

- x:

  a [character](https://rdrr.io/r/base/character.html) vector or a
  [data.frame](https://rdrr.io/r/base/data.frame.html) with one or two
  columns

- Becker:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  staphylococci should be categorised into coagulase-negative
  staphylococci ("CoNS") and coagulase-positive staphylococci ("CoPS")
  instead of their own species, according to Karsten Becker *et al.*
  (see *Source*). Please see *Details* for a full list of staphylococcal
  species that will be converted.

  This excludes *Staphylococcus aureus* at default, use `Becker = "all"`
  to also categorise *S. aureus* as "CoPS".

- Lancefield:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  beta-haemolytic *Streptococcus* should be categorised into Lancefield
  groups instead of their own species, according to Rebecca C.
  Lancefield (see *Source*). These streptococci will be categorised in
  their first group, e.g. *Streptococcus dysgalactiae* will be group C,
  although officially it was also categorised into groups G and L. .
  Please see *Details* for a full list of streptococcal species that
  will be converted.

  This excludes enterococci at default (who are in group D), use
  `Lancefield = "all"` to also categorise all enterococci as group D.

- minimum_matching_score:

  a numeric value to set as the lower limit for the [MO matching
  score](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md).
  When left blank, this will be determined automatically based on the
  character length of `x`, its [taxonomic
  kingdom](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  and [human
  pathogenicity](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md).

- keep_synonyms:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate if old,
  previously valid taxonomic names must be preserved and not be
  corrected to currently accepted names. The default is `FALSE`, which
  will return a note if old taxonomic names were processed. The default
  can be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_keep_synonyms`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md),
  i.e. `options(AMR_keep_synonyms = TRUE)` or
  `options(AMR_keep_synonyms = FALSE)`.

- reference_df:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) to be used for
  extra reference when translating `x` to a valid `mo`. See
  [`set_mo_source()`](https://amr-for-r.org/release/2.1.1/reference/mo_source.md)
  and
  [`get_mo_source()`](https://amr-for-r.org/release/2.1.1/reference/mo_source.md)
  to automate the usage of your own codes (e.g. used in your analysis or
  organisation).

- ignore_pattern:

  a Perl-compatible [regular
  expression](https://rdrr.io/r/base/regex.html) (case-insensitive) of
  which all matches in `x` must return `NA`. This can be convenient to
  exclude known non-relevant input and can also be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_ignore_pattern`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md),
  e.g.
  `options(AMR_ignore_pattern = "(not reported|contaminated flora)")`.

- cleaning_regex:

  a Perl-compatible [regular
  expression](https://rdrr.io/r/base/regex.html) (case-insensitive) to
  clean the input of `x`. Every matched part in `x` will be removed. At
  default, this is the outcome of `mo_cleaning_regex()`, which removes
  texts between brackets and texts such as "species" and "serovar". The
  default can be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_cleaning_regex`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md).

- language:

  language to translate text like "no growth", which defaults to the
  system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/2.1.1/reference/translate.md))

- info:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate if a
  progress bar should be printed if more than 25 items are to be
  coerced - the default is `TRUE` only in interactive mode

- ...:

  other arguments passed on to functions

## Value

A [character](https://rdrr.io/r/base/character.html)
[vector](https://rdrr.io/r/base/vector.html) with additional class `mo`

## Details

A microorganism (MO) code from this package (class: `mo`) is human
readable and typically looks like these examples:

      Code               Full name
      ---------------    --------------------------------------
      B_KLBSL            Klebsiella
      B_KLBSL_PNMN       Klebsiella pneumoniae
      B_KLBSL_PNMN_RHNS  Klebsiella pneumoniae rhinoscleromatis
      |   |    |    |
      |   |    |    |
      |   |    |    \---> subspecies, a 3-5 letter acronym
      |   |    \----> species, a 3-6 letter acronym
      |   \----> genus, a 4-8 letter acronym
      \----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
                                F (Fungi), PL (Plantae), P (Protozoa)

Values that cannot be coerced will be considered 'unknown' and will be
returned as the MO code `UNKNOWN` with a warning.

Use the
[`mo_*`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)
functions to get properties based on the returned code, see *Examples*.

The `as.mo()` function uses a novel [matching score
algorithm](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md)
(see *Matching Score for Microorganisms* below) to match input against
the [available microbial
taxonomy](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
in this package. This will lead to the effect that e.g. `"E. coli"` (a
microorganism highly prevalent in humans) will return the microbial ID
of *Escherichia coli* and not *Entamoeba coli* (a microorganism less
prevalent in humans), although the latter would alphabetically come
first.

With `Becker = TRUE`, the following 85 staphylococci will be converted
to the **coagulase-negative group**: *S. argensis*, *S. arlettae*, *S.
auricularis*, *S. borealis*, *S. caeli*, *S. caledonicus*, *S. canis*,
*S. capitis*, *S. capitis capitis*, *S. capitis urealyticus*, *S.
capitis ureolyticus*, *S. caprae*, *S. carnosus*, *S. carnosus
carnosus*, *S. carnosus utilis*, *S. casei*, *S. caseolyticus*, *S.
chromogenes*, *S. cohnii*, *S. cohnii cohnii*, *S. cohnii urealyticum*,
*S. cohnii urealyticus*, *S. condimenti*, *S. croceilyticus*, *S.
debuckii*, *S. devriesei*, *S. durrellii*, *S. edaphicus*, *S.
epidermidis*, *S. equorum*, *S. equorum equorum*, *S. equorum linens*,
*S. felis*, *S. fleurettii*, *S. gallinarum*, *S. haemolyticus*, *S.
hominis*, *S. hominis hominis*, *S. hominis novobiosepticus*, *S.
jettensis*, *S. kloosii*, *S. lentus*, *S. lloydii*, *S. lugdunensis*,
*S. massiliensis*, *S. microti*, *S. muscae*, *S. nepalensis*, *S.
pasteuri*, *S. petrasii*, *S. petrasii croceilyticus*, *S. petrasii
jettensis*, *S. petrasii petrasii*, *S. petrasii pragensis*, *S.
pettenkoferi*, *S. piscifermentans*, *S. pragensis*, *S. pseudoxylosus*,
*S. pulvereri*, *S. ratti*, *S. rostri*, *S. saccharolyticus*, *S.
saprophyticus*, *S. saprophyticus bovis*, *S. saprophyticus
saprophyticus*, *S. schleiferi*, *S. schleiferi schleiferi*, *S.
sciuri*, *S. sciuri carnaticus*, *S. sciuri lentus*, *S. sciuri
rodentium*, *S. sciuri sciuri*, *S. simulans*, *S. stepanovicii*, *S.
succinus*, *S. succinus casei*, *S. succinus succinus*, *S.
taiwanensis*, *S. urealyticus*, *S. ureilyticus*, *S. veratri*, *S.
vitulinus*, *S. vitulus*, *S. warneri*, and *S. xylosus*.  
The following 16 staphylococci will be converted to the
**coagulase-positive group**: *S. agnetis*, *S. argenteus*, *S.
coagulans*, *S. cornubiensis*, *S. delphini*, *S. hyicus*, *S. hyicus
chromogenes*, *S. hyicus hyicus*, *S. intermedius*, *S. lutrae*, *S.
pseudintermedius*, *S. roterodami*, *S. schleiferi coagulans*, *S.
schweitzeri*, *S. simiae*, and *S. singaporensis*.

With `Lancefield = TRUE`, the following streptococci will be converted
to their corresponding Lancefield group: *S. agalactiae* (Group B), *S.
anginosus anginosus* (Group F), *S. anginosus whileyi* (Group F), *S.
anginosus* (Group F), *S. canis* (Group G), *S. dysgalactiae
dysgalactiae* (Group C), *S. dysgalactiae equisimilis* (Group C), *S.
dysgalactiae* (Group C), *S. equi equi* (Group C), *S. equi ruminatorum*
(Group C), *S. equi zooepidemicus* (Group C), *S. equi* (Group C), *S.
pyogenes* (Group A), *S. salivarius salivarius* (Group K), *S.
salivarius thermophilus* (Group K), *S. salivarius* (Group K), and *S.
sanguinis* (Group H).

### Coping with Uncertain Results

Results of non-exact taxonomic input are based on their [matching
score](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md).
The lowest allowed score can be set with the `minimum_matching_score`
argument. At default this will be determined based on the character
length of the input, and the [taxonomic
kingdom](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
and [human
pathogenicity](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md)
of the taxonomic outcome. If values are matched with uncertainty, a
message will be shown to suggest the user to evaluate the results with
`mo_uncertainties()`, which returns a
[data.frame](https://rdrr.io/r/base/data.frame.html) with all
specifications.

To increase the quality of matching, the `cleaning_regex` argument can
be used to clean the input (i.e., `x`). This must be a [regular
expression](https://rdrr.io/r/base/regex.html) that matches parts of the
input that should be removed before the input is matched against the
[available microbial
taxonomy](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md).
It will be matched Perl-compatible and case-insensitive. The default
value of `cleaning_regex` is the outcome of the helper function
`mo_cleaning_regex()`.

There are three helper functions that can be run after using the
`as.mo()` function:

- Use `mo_uncertainties()` to get a
  [data.frame](https://rdrr.io/r/base/data.frame.html) that prints in a
  pretty format with all taxonomic names that were guessed. The output
  contains the matching score for all matches (see *Matching Score for
  Microorganisms* below).

- Use `mo_failures()` to get a
  [character](https://rdrr.io/r/base/character.html)
  [vector](https://rdrr.io/r/base/vector.html) with all values that
  could not be coerced to a valid value.

- Use `mo_renamed()` to get a
  [data.frame](https://rdrr.io/r/base/data.frame.html) with all values
  that could be coerced based on old, previously accepted taxonomic
  names.

### Microbial Prevalence of Pathogens in Humans

The coercion rules consider the prevalence of microorganisms in humans,
which is available as the `prevalence` column in the
[microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
data set. The grouping into human pathogenic prevalence is explained in
the section *Matching Score for Microorganisms* below.

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
    staphylococci.** *Expert Rev Anti Infect Ther.* 18(4):349-366;
    [doi:10.1080/14787210.2020.1730813](https://doi.org/10.1080/14787210.2020.1730813)

5.  Lancefield RC (1933). **A serological differentiation of human and
    other groups of hemolytic streptococci.** *J Exp Med.* 57(4):
    571-95;
    [doi:10.1084/jem.57.4.571](https://doi.org/10.1084/jem.57.4.571)

6.  Berends MS *et al.* (2022). **Trends in Occurrence and Phenotypic
    Resistance of Coagulase-Negative Staphylococci (CoNS) Found in Human
    Blood in the Northern Netherlands between 2013 and 2019/**
    *Micro.rganisms* 10(9), 1801;
    [doi:10.3390/microorganisms10091801](https://doi.org/10.3390/microorganisms10091801)

7.  Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing
    in Nomenclature (LPSN) moves to the DSMZ.** International Journal of
    Systematic and Evolutionary Microbiology, 70, 5607-5612;
    [doi:10.1099/ijsem.0.004332](https://doi.org/10.1099/ijsem.0.004332)
    . Accessed from <https://lpsn.dsmz.de> on December 11th, 2022.

8.  GBIF Secretariat (2022). GBIF Backbone Taxonomy. Checklist dataset
    [doi:10.15468/39omei](https://doi.org/10.15468/39omei) . Accessed
    from <https://www.gbif.org> on December 11th, 2022.

9.  Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base
    for standardized bacterial and archaeal data.** Nucleic Acids Res.,
    50(D1):D741-D74;
    [doi:10.1093/nar/gkab961](https://doi.org/10.1093/nar/gkab961) .
    Accessed from <https://bacdive.dsmz.de> on May 12th, 2023.

10. Public Health Information Network Vocabulary Access and Distribution
    System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
    Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009
    (v12). URL: <https://phinvads.cdc.gov>

11. Bartlett A *et al.* (2022). **A comprehensive list of bacterial
    pathogens infecting humans** *Microbiology* 168:001269;
    [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## Matching Score for Microorganisms

With ambiguous user input in `as.mo()` and all the
[`mo_*`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)
functions, the returned results are chosen based on their matching score
using
[`mo_matching_score()`](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md).
This matching score \\m\\, is calculated as:

![mo matching score](figures/mo_matching_score.png)

where:

- \\x\\ is the user input;

- \\n\\ is a taxonomic name (genus, species, and subspecies);

- \\l_n\\ is the length of \\n\\;

- \\lev\\ is the [Levenshtein distance
  function](https://en.wikipedia.org/wiki/Levenshtein_distance)
  (counting any insertion as 1, and any deletion or substitution as 2)
  that is needed to change \\x\\ into \\n\\;

- \\p_n\\ is the human pathogenic prevalence group of \\n\\, as
  described below;

- \\k_n\\ is the taxonomic kingdom of \\n\\, set as Bacteria = 1, Fungi
  = 1.25, Protozoa = 1.5, Archaea = 2, others = 3.

The grouping into human pathogenic prevalence \\p\\ is based on recent
work from Bartlett *et al.* (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ) who
extensively studied medical-scientific literature to categorise all
bacterial species into these groups:

- **Established**, if a taxonomic species has infected at least three
  persons in three or more references. These records have
  `prevalence = 1.0` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set;

- **Putative**, if a taxonomic species has fewer than three known cases.
  These records have `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set.

Furthermore,

- Any genus present in the **established** list also has
  `prevalence = 1.0` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set;

- Any other genus present in the **putative** list has
  `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set;

- Any other species or subspecies of which the genus is present in the
  two aforementioned groups, has `prevalence = 1.5` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set;

- Any *non-bacterial* genus, species or subspecies of which the genus is
  present in the following list, has `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set: *Absidia*, *Acanthamoeba*, *Acremonium*, *Aedes*,
  *Alternaria*, *Amoeba*, *Ancylostoma*, *Angiostrongylus*, *Anisakis*,
  *Anopheles*, *Apophysomyces*, *Aspergillus*, *Aureobasidium*,
  *Basidiobolus*, *Beauveria*, *Blastocystis*, *Blastomyces*, *Candida*,
  *Capillaria*, *Chaetomium*, *Chrysonilia*, *Cladophialophora*,
  *Cladosporium*, *Conidiobolus*, *Contracaecum*, *Cordylobia*,
  *Cryptococcus*, *Curvularia*, *Demodex*, *Dermatobia*, *Dientamoeba*,
  *Diphyllobothrium*, *Dirofilaria*, *Echinostoma*, *Entamoeba*,
  *Enterobius*, *Exophiala*, *Exserohilum*, *Fasciola*, *Fonsecaea*,
  *Fusarium*, *Giardia*, *Haloarcula*, *Halobacterium*, *Halococcus*,
  *Hendersonula*, *Heterophyes*, *Histomonas*, *Histoplasma*,
  *Hymenolepis*, *Hypomyces*, *Hysterothylacium*, *Leishmania*,
  *Malassezia*, *Malbranchea*, *Metagonimus*, *Meyerozyma*,
  *Microsporidium*, *Microsporum*, *Mortierella*, *Mucor*,
  *Mycocentrospora*, *Necator*, *Nectria*, *Ochroconis*,
  *Oesophagostomum*, *Oidiodendron*, *Opisthorchis*, *Pediculus*,
  *Penicillium*, *Phlebotomus*, *Phoma*, *Pichia*, *Piedraia*,
  *Pithomyces*, *Pityrosporum*, *Pneumocystis*, *Pseudallescheria*,
  *Pseudoterranova*, *Pulex*, *Rhizomucor*, *Rhizopus*, *Rhodotorula*,
  *Saccharomyces*, *Sarcoptes*, *Scolecobasidium*, *Scopulariopsis*,
  *Scytalidium*, *Spirometra*, *Sporobolomyces*, *Stachybotrys*,
  *Strongyloides*, *Syngamus*, *Taenia*, *Talaromyces*, *Toxocara*,
  *Trichinella*, *Trichobilharzia*, *Trichoderma*, *Trichomonas*,
  *Trichophyton*, *Trichosporon*, *Trichostrongylus*, *Trichuris*,
  *Tritirachium*, *Trombicula*, *Trypanosoma*, *Tunga*, or *Wuchereria*;

- All other records have `prevalence = 2.0` in the
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
  data set.

When calculating the matching score, all characters in \\x\\ and \\n\\
are ignored that are other than A-Z, a-z, 0-9, spaces and parentheses.

All matches are sorted descending on their matching score and for all
user input values, the top match will be returned. This will lead to the
effect that e.g., `"E. coli"` will return the microbial ID of
*Escherichia coli* (\\m = 0.688\\, a highly prevalent microorganism
found in humans) and not *Entamoeba coli* (\\m = 0.381\\, a less
prevalent microorganism in humans), although the latter would
alphabetically come first.

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

[microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)
for the [data.frame](https://rdrr.io/r/base/data.frame.html) that is
being used to determine ID's.

The
[`mo_*`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)
functions (such as
[`mo_genus()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md),
[`mo_gramstain()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md))
to get properties based on the returned code.

## Examples

``` r
# \donttest{
# These examples all return "B_STPHY_AURS", the ID of S. aureus:
as.mo(c(
  "sau", # WHONET code
  "stau",
  "STAU",
  "staaur",
  "S. aureus",
  "S aureus",
  "Sthafilokkockus aureus", # handles incorrect spelling
  "Staphylococcus aureus (MRSA)",
  "MRSA", # Methicillin Resistant S. aureus
  "VISA", # Vancomycin Intermediate S. aureus
  "VRSA", # Vancomycin Resistant S. aureus
  115329001 # SNOMED CT code
))

# Dyslexia is no problem - these all work:
as.mo(c(
  "Ureaplasma urealyticum",
  "Ureaplasma urealyticus",
  "Ureaplasmium urealytica",
  "Ureaplazma urealitycium"
))

# input will get cleaned up with the input given in the `cleaning_regex` argument,
# which defaults to `mo_cleaning_regex()`:
cat(mo_cleaning_regex(), "\n")

as.mo("Streptococcus group A")

as.mo("S. epidermidis") # will remain species: B_STPHY_EPDR
as.mo("S. epidermidis", Becker = TRUE) # will not remain species: B_STPHY_CONS

as.mo("S. pyogenes") # will remain species: B_STRPT_PYGN
as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRPA

# All mo_* functions use as.mo() internally too (see ?mo_property):
mo_genus("E. coli")
mo_gramstain("ESCO")
mo_is_intrinsic_resistant("ESCCOL", ab = "vanco")
# }
```

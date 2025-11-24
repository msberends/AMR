# Transform Arbitrary Input to Valid Microbial Taxonomy

Use this function to get a valid microorganism code (`mo`) based on
arbitrary user input. Determination is done using intelligent rules and
the complete taxonomic tree of the kingdoms Animalia, Archaea, Bacteria,
Chromista, and Protozoa, and most microbial species from the kingdom
Fungi (see *Source*). The input can be almost anything: a full name
(like `"Staphylococcus aureus"`), an abbreviated name (such as
`"S. aureus"`), an abbreviation known in the field (such as `"MRSA"`),
or just a genus. See *Examples*.

## Usage

``` r
as.mo(x, Becker = FALSE, Lancefield = FALSE,
  minimum_matching_score = NULL,
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE),
  reference_df = get_mo_source(),
  ignore_pattern = getOption("AMR_ignore_pattern", NULL),
  cleaning_regex = getOption("AMR_cleaning_regex", mo_cleaning_regex()),
  only_fungi = getOption("AMR_only_fungi", FALSE),
  language = get_AMR_locale(), info = interactive(), ...)

is.mo(x)

mo_uncertainties()

mo_renamed()

mo_failures()

mo_reset_session()

mo_cleaning_regex()
```

## Arguments

- x:

  A [character](https://rdrr.io/r/base/character.html) vector or a
  [data.frame](https://rdrr.io/r/base/data.frame.html) with one or two
  columns.

- Becker:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  staphylococci should be categorised into coagulase-negative
  staphylococci ("CoNS") and coagulase-positive staphylococci ("CoPS")
  instead of their own species, according to Karsten Becker *et al.*
  (see *Source*). Please see *Details* for a full list of staphylococcal
  species that will be converted.

  This excludes *Staphylococcus aureus* at default, use `Becker = "all"`
  to also categorise *S. aureus* as "CoPS".

- Lancefield:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
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

  A numeric value to set as the lower limit for the [MO matching
  score](https://amr-for-r.org/reference/mo_matching_score.md). When
  left blank, this will be determined automatically based on the
  character length of `x`, its [taxonomic
  kingdom](https://amr-for-r.org/reference/microorganisms.md) and [human
  pathogenicity](https://amr-for-r.org/reference/mo_matching_score.md).

- keep_synonyms:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if old,
  previously valid taxonomic names must be preserved and not be
  corrected to currently accepted names. The default is `FALSE`, which
  will return a note if old taxonomic names were processed. The default
  can be set with the package option
  [`AMR_keep_synonyms`](https://amr-for-r.org/reference/AMR-options.md),
  i.e. `options(AMR_keep_synonyms = TRUE)` or
  `options(AMR_keep_synonyms = FALSE)`.

- reference_df:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) to be used for
  extra reference when translating `x` to a valid `mo`. See
  [`set_mo_source()`](https://amr-for-r.org/reference/mo_source.md) and
  [`get_mo_source()`](https://amr-for-r.org/reference/mo_source.md) to
  automate the usage of your own codes (e.g. used in your analysis or
  organisation).

- ignore_pattern:

  A Perl-compatible [regular
  expression](https://rdrr.io/r/base/regex.html) (case-insensitive) of
  which all matches in `x` must return `NA`. This can be convenient to
  exclude known non-relevant input and can also be set with the package
  option
  [`AMR_ignore_pattern`](https://amr-for-r.org/reference/AMR-options.md),
  e.g.
  `options(AMR_ignore_pattern = "(not reported|contaminated flora)")`.

- cleaning_regex:

  A Perl-compatible [regular
  expression](https://rdrr.io/r/base/regex.html) (case-insensitive) to
  clean the input of `x`. Every matched part in `x` will be removed. At
  default, this is the outcome of `mo_cleaning_regex()`, which removes
  texts between brackets and texts such as "species" and "serovar". The
  default can be set with the package option
  [`AMR_cleaning_regex`](https://amr-for-r.org/reference/AMR-options.md).

- only_fungi:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if only
  fungi must be found, making sure that e.g. misspellings always return
  records from the kingdom of Fungi. This can be set globally for [all
  microorganism
  functions](https://amr-for-r.org/reference/mo_property.md) with the
  package option
  [`AMR_only_fungi`](https://amr-for-r.org/reference/AMR-options.md),
  i.e. `options(AMR_only_fungi = TRUE)`.

- language:

  Language to translate text like "no growth", which defaults to the
  system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)).

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate that info
  must be printed, e.g. a progress bar when more than 25 items are to be
  coerced, or a list with old taxonomic names. The default is `TRUE`
  only in interactive mode.

- ...:

  Other arguments passed on to functions.

## Value

A [character](https://rdrr.io/r/base/character.html)
[vector](https://rdrr.io/r/base/vector.html) with additional class `mo`

## Details

A microorganism (MO) code from this package (class: `mo`) is
human-readable and typically looks like these examples:

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
      \----> kingdom: A (Archaea), AN (Animalia), B (Bacteria),
                      C (Chromista), F (Fungi), PL (Plantae),
                      P (Protozoa)

Values that cannot be coerced will be considered 'unknown' and will
return the MO code `UNKNOWN` with a warning.

Use the [`mo_*`](https://amr-for-r.org/reference/mo_property.md)
functions to get properties based on the returned code, see *Examples*.

The `as.mo()` function uses a novel and scientifically validated
([doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03) )
matching score algorithm (see *Matching Score for Microorganisms* below)
to match input against the [available microbial
taxonomy](https://amr-for-r.org/reference/microorganisms.md) in this
package. This implicates that e.g. `"E. coli"` (a microorganism highly
prevalent in humans) will return the microbial ID of *Escherichia coli*
and not *Entamoeba coli* (a microorganism less prevalent in humans),
although the latter would alphabetically come first.

### Coping with Uncertain Results

Results of non-exact taxonomic input are based on their [matching
score](https://amr-for-r.org/reference/mo_matching_score.md). The lowest
allowed score can be set with the `minimum_matching_score` argument. At
default this will be determined based on the character length of the
input, the [taxonomic
kingdom](https://amr-for-r.org/reference/microorganisms.md), and the
[human
pathogenicity](https://amr-for-r.org/reference/mo_matching_score.md) of
the taxonomic outcome. If values are matched with uncertainty, a message
will be shown to suggest the user to inspect the results with
`mo_uncertainties()`, which returns a
[data.frame](https://rdrr.io/r/base/data.frame.html) with all
specifications.

To increase the quality of matching, the `cleaning_regex` argument is
used to clean the input. This must be a [regular
expression](https://rdrr.io/r/base/regex.html) that matches parts of the
input that should be removed before the input is matched against the
[available microbial
taxonomy](https://amr-for-r.org/reference/microorganisms.md). It will be
matched Perl-compatible and case-insensitive. The default value of
`cleaning_regex` is the outcome of the helper function
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

### For Mycologists

The [matching score
algorithm](https://amr-for-r.org/reference/mo_matching_score.md) gives
precedence to bacteria over fungi. If you are only analysing fungi, be
sure to use `only_fungi = TRUE`, or better yet, add this to your code
and run it once every session:

    options(AMR_only_fungi = TRUE)

This will make sure that no bacteria or other 'non-fungi' will be
returned by `as.mo()`, or any of the
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions.

### Coagulase-negative and Coagulase-positive Staphylococci

With `Becker = TRUE`, the following staphylococci will be converted to
their corresponding coagulase group:

- Coagulase-negative: *S. americanisciuri*, *S. argensis*, *S.
  arlettae*, *S. auricularis*, *S. borealis*, *S. brunensis*, *S.
  caeli*, *S. caledonicus*, *S. canis*, *S. capitis*, *S. capitis
  capitis*, *S. capitis urealyticus*, *S. capitis ureolyticus*, *S.
  caprae*, *S. carnosus*, *S. carnosus carnosus*, *S. carnosus utilis*,
  *S. casei*, *S. caseolyticus*, *S. chromogenes*, *S. cohnii*, *S.
  cohnii cohnii*, *S. cohnii urealyticum*, *S. cohnii urealyticus*, *S.
  condimenti*, *S. croceilyticus*, *S. debuckii*, *S. devriesei*, *S.
  durrellii*, *S. edaphicus*, *S. epidermidis*, *S. equorum*, *S.
  equorum equorum*, *S. equorum linens*, *S. felis*, *S. fleurettii*,
  *S. gallinarum*, *S. haemolyticus*, *S. hominis*, *S. hominis
  hominis*, *S. hominis novobiosepticus*, *S. jettensis*, *S. kloosii*,
  *S. lentus*, *S. lloydii*, *S. lugdunensis*, *S. marylandisciuri*, *S.
  massiliensis*, *S. microti*, *S. muscae*, *S. nepalensis*, *S.
  pasteuri*, *S. petrasii*, *S. petrasii croceilyticus*, *S. petrasii
  jettensis*, *S. petrasii petrasii*, *S. petrasii pragensis*, *S.
  pettenkoferi*, *S. piscifermentans*, *S. pragensis*, *S.
  pseudoxylosus*, *S. pulvereri*, *S. ratti*, *S. rostri*, *S.
  saccharolyticus*, *S. saprophyticus*, *S. saprophyticus bovis*, *S.
  saprophyticus saprophyticus*, *S. schleiferi*, *S. schleiferi
  schleiferi*, *S. sciuri*, *S. sciuri carnaticus*, *S. sciuri lentus*,
  *S. sciuri rodentium*, *S. sciuri sciuri*, *S. shinii*, *S. simulans*,
  *S. stepanovicii*, *S. succinus*, *S. succinus casei*, *S. succinus
  succinus*, *S. taiwanensis*, *S. urealyticus*, *S. ureilyticus*, *S.
  veratri*, *S. vitulinus*, *S. vitulus*, *S. warneri*, and *S. xylosus*

- Coagulase-positive: *S. agnetis*, *S. argenteus*, *S. coagulans*, *S.
  cornubiensis*, *S. delphini*, *S. hyicus*, *S. hyicus chromogenes*,
  *S. hyicus hyicus*, *S. intermedius*, *S. lutrae*, *S.
  pseudintermedius*, *S. roterodami*, *S. schleiferi coagulans*, *S.
  schweitzeri*, *S. simiae*, and *S. singaporensis*

This is based on:

- Becker K *et al.* (2014). **Coagulase-Negative Staphylococci.** *Clin
  Microbiol Rev.* 27(4): 870-926;
  [doi:10.1128/CMR.00109-13](https://doi.org/10.1128/CMR.00109-13)

- Becker K *et al.* (2019). **Implications of identifying the recently
  defined members of the *S. aureus* complex, *S. argenteus* and *S.
  schweitzeri*: A position paper of members of the ESCMID Study Group
  for staphylococci and Staphylococcal Diseases (ESGS).** *Clin
  Microbiol Infect*;
  [doi:10.1016/j.cmi.2019.02.028](https://doi.org/10.1016/j.cmi.2019.02.028)

- Becker K *et al.* (2020). **Emergence of coagulase-negative
  staphylococci.** *Expert Rev Anti Infect Ther.* 18(4):349-366;
  [doi:10.1080/14787210.2020.1730813](https://doi.org/10.1080/14787210.2020.1730813)

For newly named staphylococcal species, such as *S. brunensis* (2024)
and *S. shinii* (2023), we looked up the scientific reference to make
sure the species are considered for the correct coagulase group.

### Lancefield Groups in Streptococci

With `Lancefield = TRUE`, the following streptococci will be converted
to their corresponding Lancefield group:

- Streptococcus Group A: *S. pyogenes*

- Streptococcus Group B: *S. agalactiae*

- Streptococcus Group C: *S. dysgalactiae*, *S. dysgalactiae
  dysgalactiae*, *S. dysgalactiae equisimilis*, *S. equi*, *S. equi
  equi*, *S. equi ruminatorum*, and *S. equi zooepidemicus*

- Streptococcus Group F: *S. anginosus*, *S. anginosus anginosus*, *S.
  anginosus whileyi*, *S. constellatus*, *S. constellatus constellatus*,
  *S. constellatus pharyngis*, *S. constellatus viborgensis*, and *S.
  intermedius*

- Streptococcus Group G: *S. canis*, *S. dysgalactiae*, *S. dysgalactiae
  dysgalactiae*, and *S. dysgalactiae equisimilis*

- Streptococcus Group H: *S. sanguinis*

- Streptococcus Group K: *S. salivarius*, *S. salivarius salivarius*,
  and *S. salivarius thermophilus*

- Streptococcus Group L: *S. dysgalactiae*, *S. dysgalactiae
  dysgalactiae*, and *S. dysgalactiae equisimilis*

This is based on:

- Lancefield RC (1933). **A serological differentiation of human and
  other groups of hemolytic streptococci.** *J Exp Med.* 57(4): 571-95;
  [doi:10.1084/jem.57.4.571](https://doi.org/10.1084/jem.57.4.571)

## Source

- Berends MS *et al.* (2022). **AMR: An R Package for Working with
  Antimicrobial Resistance Data**. *Journal of Statistical Software*,
  104(3), 1-31;
  [doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)

- Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing
  in Nomenclature (LPSN) moves to the DSMZ.** International Journal of
  Systematic and Evolutionary Microbiology, 70, 5607-5612;
  [doi:10.1099/ijsem.0.004332](https://doi.org/10.1099/ijsem.0.004332) .
  Accessed from <https://lpsn.dsmz.de> on June 24th, 2024.

- Vincent, R *et al* (2013). **MycoBank gearing up for new horizons.**
  IMA Fungus, 4(2), 371-9;
  [doi:10.5598/imafungus.2013.04.02.16](https://doi.org/10.5598/imafungus.2013.04.02.16)
  . Accessed from <https://www.mycobank.org> on June 24th, 2024.

- GBIF Secretariat (2023). GBIF Backbone Taxonomy. Checklist dataset
  [doi:10.15468/39omei](https://doi.org/10.15468/39omei) . Accessed from
  <https://www.gbif.org> on June 24th, 2024.

- Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base
  for standardized bacterial and archaeal data.** Nucleic Acids Res.,
  50(D1):D741-D74;
  [doi:10.1093/nar/gkab961](https://doi.org/10.1093/nar/gkab961) .
  Accessed from <https://bacdive.dsmz.de> on July 16th, 2024.

- Public Health Information Network Vocabulary Access and Distribution
  System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
  Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).
  URL: <https://www.cdc.gov/phin/php/phinvads/>

- Bartlett A *et al.* (2022). **A comprehensive list of bacterial
  pathogens infecting humans** *Microbiology* 168:001269;
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## Matching Score for Microorganisms

With ambiguous user input in `as.mo()` and all the
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions, the
returned results are chosen based on their matching score using
[`mo_matching_score()`](https://amr-for-r.org/reference/mo_matching_score.md).
This matching score \\m\\, is calculated as:

\$\$m\_{(x, n)} = \frac{l\_{n} - 0.5 \cdot \min \begin{cases}l\_{n} \\
\textrm{lev}(x, n)\end{cases}}{l\_{n} \cdot p\_{n} \cdot k\_{n}}\$\$

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
  = 1.25, Protozoa = 1.5, Chromista = 1.75, Archaea = 2, others = 3.

The grouping into human pathogenic prevalence \\p\\ is based on recent
work from Bartlett *et al.* (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ) who
extensively studied medical-scientific literature to categorise all
bacterial species into these groups:

- **Established**, if a taxonomic species has infected at least three
  persons in three or more references. These records have
  `prevalence = 1.15` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set;

- **Putative**, if a taxonomic species has fewer than three known cases.
  These records have `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set.

Furthermore,

- Genera from the World Health Organization's (WHO) Priority Pathogen
  List have `prevalence = 1.0` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set;

- Any genus present in the **established** list also has
  `prevalence = 1.15` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set;

- Any other genus present in the **putative** list has
  `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set;

- Any other species or subspecies of which the genus is present in the
  two aforementioned groups, has `prevalence = 1.5` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set;

- Any *non-bacterial* genus, species or subspecies of which the genus is
  present in the following list, has `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set: *Absidia*, *Acanthamoeba*, *Acremonium*, *Actinomucor*,
  *Aedes*, *Alternaria*, *Amoeba*, *Ancylostoma*, *Angiostrongylus*,
  *Anisakis*, *Anopheles*, *Apophysomyces*, *Arthroderma*,
  *Aspergillus*, *Aureobasidium*, *Basidiobolus*, *Beauveria*,
  *Bipolaris*, *Blastobotrys*, *Blastocystis*, *Blastomyces*, *Candida*,
  *Capillaria*, *Chaetomium*, *Chilomastix*, *Chrysonilia*,
  *Chrysosporium*, *Cladophialophora*, *Cladosporium*, *Clavispora*,
  *Coccidioides*, *Cokeromyces*, *Conidiobolus*, *Coniochaeta*,
  *Contracaecum*, *Cordylobia*, *Cryptococcus*, *Cryptosporidium*,
  *Cunninghamella*, *Curvularia*, *Cyberlindnera*, *Debaryozyma*,
  *Demodex*, *Dermatobia*, *Dientamoeba*, *Diphyllobothrium*,
  *Dirofilaria*, *Echinostoma*, *Entamoeba*, *Enterobius*,
  *Epidermophyton*, *Exidia*, *Exophiala*, *Exserohilum*, *Fasciola*,
  *Fonsecaea*, *Fusarium*, *Geotrichum*, *Giardia*, *Graphium*,
  *Haloarcula*, *Halobacterium*, *Halococcus*, *Hansenula*,
  *Hendersonula*, *Heterophyes*, *Histomonas*, *Histoplasma*, *Hortaea*,
  *Hymenolepis*, *Hypomyces*, *Hysterothylacium*, *Kloeckera*,
  *Kluyveromyces*, *Kodamaea*, *Lacazia*, *Leishmania*, *Lichtheimia*,
  *Lodderomyces*, *Lomentospora*, *Madurella*, *Malassezia*,
  *Malbranchea*, *Metagonimus*, *Meyerozyma*, *Microsporidium*,
  *Microsporum*, *Millerozyma*, *Mortierella*, *Mucor*,
  *Mycocentrospora*, *Nannizzia*, *Necator*, *Nectria*, *Ochroconis*,
  *Oesophagostomum*, *Oidiodendron*, *Opisthorchis*, *Paecilomyces*,
  *Paracoccidioides*, *Pediculus*, *Penicillium*, *Phaeoacremonium*,
  *Phaeomoniella*, *Phialophora*, *Phlebotomus*, *Phoma*, *Pichia*,
  *Piedraia*, *Pithomyces*, *Pityrosporum*, *Pneumocystis*,
  *Pseudallescheria*, *Pseudoscopulariopsis*, *Pseudoterranova*,
  *Pulex*, *Purpureocillium*, *Quambalaria*, *Rhinocladiella*,
  *Rhizomucor*, *Rhizopus*, *Rhodotorula*, *Saccharomyces*, *Saksenaea*,
  *Saprochaete*, *Sarcoptes*, *Scedosporium*, *Schistosoma*,
  *Schizosaccharomyces*, *Scolecobasidium*, *Scopulariopsis*,
  *Scytalidium*, *Spirometra*, *Sporobolomyces*, *Sporopachydermia*,
  *Sporothrix*, *Sporotrichum*, *Stachybotrys*, *Strongyloides*,
  *Syncephalastrum*, *Syngamus*, *Taenia*, *Talaromyces*, *Teleomorph*,
  *Toxocara*, *Trichinella*, *Trichobilharzia*, *Trichoderma*,
  *Trichomonas*, *Trichophyton*, *Trichosporon*, *Trichostrongylus*,
  *Trichuris*, *Tritirachium*, *Trombicula*, *Trypanosoma*, *Tunga*,
  *Ulocladium*, *Ustilago*, *Verticillium*, *Wallemia*, *Wangiella*,
  *Wickerhamomyces*, *Wuchereria*, *Yarrowia*, or *Zygosaccharomyces*;

- All other records have `prevalence = 2.0` in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
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

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## See also

[microorganisms](https://amr-for-r.org/reference/microorganisms.md) for
the [data.frame](https://rdrr.io/r/base/data.frame.html) that is being
used to determine ID's.

The [`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions
(such as [`mo_genus()`](https://amr-for-r.org/reference/mo_property.md),
[`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md)) to
get properties based on the returned code.

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
#> Class 'mo'
#>  [1] B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS
#>  [6] B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS B_STPHY_AURS
#> [11] B_STPHY_AURS B_STPHY_AURS

# Dyslexia is no problem - these all work:
as.mo(c(
  "Ureaplasma urealyticum",
  "Ureaplasma urealyticus",
  "Ureaplasmium urealytica",
  "Ureaplazma urealitycium"
))
#> Class 'mo'
#> [1] B_URPLS_URLY B_URPLS_URLY B_URPLS_URLY B_URPLS_URLY

# input will get cleaned up with the input given in the `cleaning_regex` argument,
# which defaults to `mo_cleaning_regex()`:
cat(mo_cleaning_regex(), "\n")
#> ([^A-Za-z- \(\)\[\]{}]+|([({]|\[).+([})]|\])|(^| )( ?[a-z-]+[-](resistant|susceptible) ?|e?spp([^a-z]+|$)|e?ssp([^a-z]+|$)|serogr.?up[a-z]*|e?ss([^a-z]+|$)|e?sp([^a-z]+|$)|var([^a-z]+|$)|serovar[a-z]*|sube?species|biovar[a-z]*|e?species|Ig[ADEGM]|e?subsp|biotype|titer|dummy)) 

as.mo("Streptococcus group A")
#> Class 'mo'
#> [1] B_STRPT_GRPA

as.mo("S. epidermidis") # will remain species: B_STPHY_EPDR
#> Class 'mo'
#> [1] B_STPHY_EPDR
as.mo("S. epidermidis", Becker = TRUE) # will not remain species: B_STPHY_CONS
#> Class 'mo'
#> [1] B_STPHY_CONS

as.mo("S. pyogenes") # will remain species: B_STRPT_PYGN
#> Class 'mo'
#> [1] B_STRPT_PYGN
as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRPA
#> Class 'mo'
#> [1] B_STRPT_GRPA

# All mo_* functions use as.mo() internally too (see ?mo_property):
mo_genus("E. coli")
#> [1] "Escherichia"
mo_gramstain("ESCO")
#> [1] "Gram-negative"
mo_is_intrinsic_resistant("ESCCOL", ab = "vanco")
#> ℹ Determining intrinsic resistance based on 'EUCAST Expected Resistant
#>   Phenotypes' v1.2 (2023). This note will be shown once per session.
#> [1] TRUE
# }
```

# Data Set with 52 171 Microorganisms

A data set containing the full microbial taxonomy (**last updated:
December 11th, 2022**) of five kingdoms from the List of Prokaryotic
names with Standing in Nomenclature (LPSN) and the Global Biodiversity
Information Facility (GBIF). This data set is the backbone of this `AMR`
package. MO codes can be looked up using
[`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md).

## Usage

``` r
microorganisms
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 52
171 observations and 23 variables:

- `mo`  
  ID of microorganism as used by this package. *This is a unique
  identifier.*

- `fullname`  
  Full name, like `"Escherichia coli"`. For the taxonomic ranks genus,
  species and subspecies, this is the 'pasted' text of genus, species,
  and subspecies. For all taxonomic ranks higher than genus, this is the
  name of the taxon. *This is a unique identifier.*

- `status`  
  Status of the taxon, either "accepted" or "synonym"

- `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`,
  `subspecies`  
  Taxonomic rank of the microorganism

- `rank`  
  Text of the taxonomic rank of the microorganism, such as `"species"`
  or `"genus"`

- `ref`  
  Author(s) and year of related scientific publication. This contains
  only the *first surname* and year of the *latest* authors, e.g.
  "Wallis *et al.* 2006 *emend.* Smith and Jones 2018" becomes "Smith
  *et al.*, 2018". This field is directly retrieved from the source
  specified in the column `source`. Moreover, accents were removed to
  comply with CRAN that only allows ASCII characters.

- `lpsn`  
  Identifier ('Record number') of the List of Prokaryotic names with
  Standing in Nomenclature (LPSN). This will be the first/highest LPSN
  identifier to keep one identifier per row. For example, *Acetobacter
  ascendens* has LPSN Record number 7864 and 11011. Only the first is
  available in the `microorganisms` data set.

- `oxygen_tolerance`  
  Oxygen tolerance, either "aerobe", "anaerobe",
  "anaerobe/microaerophile", "facultative anaerobe", "likely facultative
  anaerobe", or "microaerophile". These data were retrieved from BacDive
  (see *Source*). Items that contain "likely" are missing from BacDive
  and were extrapolated from other species within the same genus to
  guess the oxygen tolerance. Currently 73.4% of all ~37 000 bacteria in
  the data set contain an oxygen tolerance.

- `lpsn_parent`  
  LPSN identifier of the parent taxon

- `lpsn_renamed_to`  
  LPSN identifier of the currently valid taxon

- `gbif`  
  Identifier ('taxonID') of the Global Biodiversity Information Facility
  (GBIF)

- `gbif_parent`  
  GBIF identifier of the parent taxon

- `gbif_renamed_to`  
  GBIF identifier of the currently valid taxon

- `source`  
  Either "GBIF", "LPSN", or "manually added" (see *Source*)

- `prevalence`  
  Prevalence of the microorganism according to Bartlett *et al.* (2022,
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ),
  see
  [`mo_matching_score()`](https://amr-for-r.org/release/2.1.1/reference/mo_matching_score.md)
  for the full explanation

- `snomed`  
  Systematized Nomenclature of Medicine (SNOMED) code of the
  microorganism, version of July 1st, 2021 (see *Source*). Use
  [`mo_snomed()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)
  to retrieve it quickly, see
  [`mo_property()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md).

## Source

- Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing
  in Nomenclature (LPSN) moves to the DSMZ.** International Journal of
  Systematic and Evolutionary Microbiology, 70, 5607-5612;
  [doi:10.1099/ijsem.0.004332](https://doi.org/10.1099/ijsem.0.004332) .
  Accessed from <https://lpsn.dsmz.de> on December 11th, 2022.

- GBIF Secretariat (2022). GBIF Backbone Taxonomy. Checklist dataset
  [doi:10.15468/39omei](https://doi.org/10.15468/39omei) . Accessed from
  <https://www.gbif.org> on December 11th, 2022.

- Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base
  for standardized bacterial and archaeal data.** Nucleic Acids Res.,
  50(D1):D741-D74;
  [doi:10.1093/nar/gkab961](https://doi.org/10.1093/nar/gkab961) .
  Accessed from <https://bacdive.dsmz.de> on May 12th, 2023.

- Public Health Information Network Vocabulary Access and Distribution
  System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
  Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).
  URL: <https://phinvads.cdc.gov>

- Grimont *et al.* (2007). Antigenic Formulae of the Salmonella
  Serovars, 9th Edition. WHO Collaborating Centre for Reference and
  Research on *Salmonella* (WHOCC-SALM).

- Bartlett *et al.* (2022). **A comprehensive list of bacterial
  pathogens infecting humans** *Microbiology* 168:001269;
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## Details

Please note that entries are only based on the List of Prokaryotic names
with Standing in Nomenclature (LPSN) and the Global Biodiversity
Information Facility (GBIF) (see below). Since these sources incorporate
entries based on (recent) publications in the International Journal of
Systematic and Evolutionary Microbiology (IJSEM), it can happen that the
year of publication is sometimes later than one might expect.

For example, *Staphylococcus pettenkoferi* was described for the first
time in Diagnostic Microbiology and Infectious Disease in 2002
([doi:10.1016/s0732-8893(02)00399-1](https://doi.org/10.1016/s0732-8893%2802%2900399-1)
), but it was not before 2007 that a publication in IJSEM followed
([doi:10.1099/ijs.0.64381-0](https://doi.org/10.1099/ijs.0.64381-0) ).
Consequently, the `AMR` package returns 2007 for
`mo_year("S. pettenkoferi")`.

## Included Taxa

Included taxonomic data are:

- All ~37 000 (sub)species from the kingdoms of Archaea and Bacteria

- ~7 900 (sub)species from the kingdom of Fungi. The kingdom of Fungi is
  a very large taxon with almost 300,000 different (sub)species, of
  which most are not microbial (but rather macroscopic, like mushrooms).
  Because of this, not all fungi fit the scope of this package. Only
  relevant fungi are covered (such as all species of *Aspergillus*,
  *Candida*, *Cryptococcus*, *Histoplasma*, *Pneumocystis*,
  *Saccharomyces* and *Trichophyton*).

- ~5 100 (sub)species from the kingdom of Protozoa

- ~1 400 (sub)species from 43 other relevant genera from the kingdom of
  Animalia (such as *Strongyloides* and *Taenia*)

- All ~9 800 previously accepted names of all included (sub)species
  (these were taxonomically renamed)

- The complete taxonomic tree of all included (sub)species: from kingdom
  to subspecies

- The identifier of the parent taxons

- The year and first author of the related scientific publication

### Manual additions

For convenience, some entries were added manually:

- ~1 500 entries of *Salmonella*, such as the city-like serovars and
  groups A to H

- 36 species groups (such as the beta-haemolytic *Streptococcus* groups
  A to K, coagulase-negative *Staphylococcus* (CoNS), *Mycobacterium
  tuberculosis* complex, etc.), of which the group compositions are
  stored in the
  [microorganisms.groups](https://amr-for-r.org/release/2.1.1/reference/microorganisms.groups.md)
  data set

- 1 entry of *Blastocystis* (*B. hominis*), although it officially does
  not exist (Noel *et al.* 2005, PMID 15634993)

- 1 entry of *Moraxella* (*M. catarrhalis*), which was formally named
  *Branhamella catarrhalis* (Catlin, 1970) though this change was never
  accepted within the field of clinical microbiology

- 8 other 'undefined' entries (unknown, unknown Gram-negatives, unknown
  Gram-positives, unknown yeast, unknown fungus, and unknown anaerobic
  Gram-pos/Gram-neg bacteria)

The syntax used to transform the original data to a cleansed R format,
can be found here:
<https://github.com/msberends/AMR/blob/main/data-raw/reproduction_of_microorganisms.R>.

### Direct download

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## About the Records from LPSN (see *Source*)

LPSN is the main source for bacteriological taxonomy of this `AMR`
package.

The List of Prokaryotic names with Standing in Nomenclature (LPSN)
provides comprehensive information on the nomenclature of prokaryotes.
LPSN is a free to use service founded by Jean P. Euzeby in 1997 and
later on maintained by Aidan C. Parte.

## See also

[`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md),
[`mo_property()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md),
[microorganisms.groups](https://amr-for-r.org/release/2.1.1/reference/microorganisms.groups.md),
[microorganisms.codes](https://amr-for-r.org/release/2.1.1/reference/microorganisms.codes.md),
[intrinsic_resistant](https://amr-for-r.org/release/2.1.1/reference/intrinsic_resistant.md)

## Examples

``` r
microorganisms
```

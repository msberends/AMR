# Data Set with 96 982 Taxonomic Records of Microorganisms

A data set containing the full microbial taxonomy (**last updated: 7th
of May, 2026**) of 15 kingdoms. This data set is the backbone of this
`AMR` package. MO codes can be looked up using
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and microorganism
properties can be looked up using any of the
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions.

This data set is carefully crafted, yet made 100% reproducible from
public and authoritative taxonomic sources (using [this
script](https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R)),
namely: *List of Prokaryotic names with Standing in Nomenclature (LPSN)*
for bacteria, *MycoBank* for fungi, and *Global Biodiversity Information
Facility (GBIF), via Catalogue of Life (COL)* for all others taxons.

## Usage

``` r
microorganisms
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 96
982 observations and 28 variables:

- `mo`  
  ID of microorganism as used by this package. ***This is a unique
  identifier.***

- `fullname`  
  Full name, like `"Escherichia coli"`. For the taxonomic ranks genus,
  species and subspecies, this is the 'pasted' text of genus, species,
  and subspecies. For all taxonomic ranks higher than genus, this is the
  name of the taxon. ***This is a unique identifier.***

- `status`  
  Status of the taxon, either `"accepted"`, `"synonym"`, or `"unknown"`

- `domain`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`,
  `species`, `subspecies`  
  Taxonomic rank of the microorganism. Note that for fungi, *phylum* is
  used for their taxonomic *division*. Also, for fungi, *subkingdom* and
  *subdivision* were left out since they do not occur in the bacterial
  taxonomy. For all species outside the domains of Bacteria and Archaea,
  the `domain` and `kingdom` are identical.

- `rank`  
  Text of the taxonomic rank of the microorganism, such as `"species"`
  or `"genus"`

- `ref`  
  Abbreviated authority citation for the nomenclatural act that
  established the current name combination, following ICNP conventions.
  For species described in their current genus (*sp. nov.*), this is the
  original description author(s) and year. For species transferred to a
  different genus (*comb. nov.*), this is the reclassification author(s)
  and year. Emendations are excluded. For synonyms, this is the
  authority under which the synonym was originally published. This field
  is directly retrieved from the source specified in the column
  `source`. Diacritics were removed to comply with CRAN, that only
  allows ASCII characters.

- `oxygen_tolerance`  
  Oxygen tolerance, either `"aerobe"`, `"anaerobe"`,
  `"anaerobe/microaerophile"`, `"facultative anaerobe"`,
  `"likely facultative anaerobe"`, `"microaerophile"`, or NA. These data
  were retrieved from BacDive (see *Source*). Items that contain
  "likely" are missing from BacDive and were extrapolated from other
  species within the same genus to guess the oxygen tolerance. Currently
  1.3784 × 10⁶% of all 2 bacteria in the data set contain an oxygen
  tolerance.

- `morphology`  
  Morphology (cell shape), either `"cocci"`, `"coccobacilli"`,
  `"filamentous"`, `"likely cocci"`, `"likely coccobacilli"`,
  `"likely filamentous"`, `"likely rods"`, `"likely spirilla"`,
  `"rods"`, `"spirilla"`, or NA. These data were retrieved from BacDive
  (see *Source*). Genera that are clinically established as coccobacilli
  (the HACEK group and beyond, such as *Haemophilus* and
  *Acinetobacter*) are classified as such regardless of BacDive majority
  vote. Items that contain "likely" are missing from BacDive and were
  extrapolated from other species within the same genus. Currently
  1.3232 × 10⁶% of all 2 bacteria in the data set contain a morphology.

- `source`  
  Either `"GBIF"`, `"LPSN"`, `"MycoBank"`, or `"manually added"` (see
  *Source*)

- `lpsn`  
  Identifier ('Record number') of List of Prokaryotic names with
  Standing in Nomenclature (LPSN). This will be the first/highest LPSN
  identifier to keep one identifier per row. For example, *Acetobacter
  ascendens* has LPSN Record number 7864 and 11011. Only the first is
  available in the `microorganisms` data set. ***This is a unique
  identifier***, though available for only ~36 000 records.

- `lpsn_parent`  
  LPSN identifier of the parent taxon

- `lpsn_renamed_to`  
  LPSN identifier of the currently valid taxon

- `mycobank`  
  Identifier ('MycoBank \#') of MycoBank. ***This is a unique
  identifier***, though available for only ~25 000 records.

- `mycobank_parent`  
  MycoBank identifier of the parent taxon

- `mycobank_renamed_to`  
  MycoBank identifier of the currently valid taxon

- `gbif`  
  Identifier ('taxonID') of Global Biodiversity Information Facility
  (GBIF), via Catalogue of Life (COL). ***This is a unique
  identifier***, though available for only ~79 000 records.

- `gbif_parent`  
  GBIF identifier of the parent taxon

- `gbif_renamed_to`  
  GBIF identifier of the currently valid taxon

- `prevalence`  
  Prevalence of the microorganism based on Bartlett *et al.* (2022,
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ),
  see
  [`mo_matching_score()`](https://amr-for-r.org/reference/mo_matching_score.md)
  for the full explanation

- `snomed`  
  Systematized Nomenclature of Medicine (SNOMED) code of the
  microorganism, version of 16th of July, 2024 (see *Source*). Use
  [`mo_snomed()`](https://amr-for-r.org/reference/mo_property.md) to
  retrieve it quickly, see
  [`mo_property()`](https://amr-for-r.org/reference/mo_property.md).

## Details

Please note that entries are only based on LPSN, MycoBank, and GBIF (see
below). Since these sources incorporate entries based on (recent)
publications in the International Journal of Systematic and Evolutionary
Microbiology (IJSEM), it can happen that the year of publication is
sometimes later than one might expect.

For example, *Staphylococcus pettenkoferi* was described for the first
time in Diagnostic Microbiology and Infectious Disease in 2002
([doi:10.1016/s0732-8893(02)00399-1](https://doi.org/10.1016/s0732-8893%2802%2900399-1)
), but it was not until 2007 that a publication in IJSEM followed
([doi:10.1099/ijs.0.64381-0](https://doi.org/10.1099/ijs.0.64381-0) ).
Consequently, the `AMR` package returns 2007 for
`mo_year("S. pettenkoferi")`.

## Included Taxa

Included taxonomic data from [LPSN](https://lpsn.dsmz.de),
[MycoBank](https://www.mycobank.org), and [GBIF](https://www.gbif.org)
are:

- All 2 (sub)species from the kingdoms of Archaea and Bacteria

- ~36 000 species from the kingdom of Fungi. The kingdom of Fungi is a
  very large taxon with almost 300,000 different (sub)species, of which
  most are not microbial (but rather macroscopic, like mushrooms).
  Because of this, not all fungi fit the scope of this package. Only
  relevant fungi are covered (such as all species of *Aspergillus*,
  *Candida*, *Cryptococcus*, *Histoplasma*, *Pneumocystis*,
  *Saccharomyces* and *Trichophyton*).

- ~11 000 (sub)species from the kingdom of Protozoa

- ~2 000 (sub)species from ~60 other relevant genera from the kingdom of
  Animalia (such as *Strongyloides* and *Taenia*)

- All ~31 000 previously accepted names of all included (sub)species
  (these were taxonomically renamed)

- The complete taxonomic tree of all included (sub)species: from kingdom
  to subspecies

- The identifier of the parent taxons

- The year and first author of the related scientific publication

### Manual additions

For convenience, some entries were added manually:

- ~1 500 entries of *Salmonella*, such as the city-like serovars and
  groups A to H

- 37 species groups (such as the beta-haemolytic *Streptococcus* groups
  A to K, coagulase-negative *Staphylococcus* (CoNS), *Mycobacterium
  tuberculosis* complex, etc.), of which the group compositions are
  stored in the
  [microorganisms.groups](https://amr-for-r.org/reference/microorganisms.groups.md)
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
can be [found
here](https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R).

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

## References

Taxonomic entries were imported in this order of importance:

1.  List of Prokaryotic names with Standing in Nomenclature (LPSN):  
      
    Freese, HM *et al.* (2026). **TYGS and LPSN in 2025: a Global Core
    Biodata Resource for genome-based classification and nomenclature of
    prokaryotes within DSMZ Digital Diversity.** Nucleic Acids Research,
    54, D884–D891;
    [doi:10.1093/nar/gkaf1110](https://doi.org/10.1093/nar/gkaf1110) .
    Accessed from <https://lpsn.dsmz.de> on 7th of May, 2026.

2.  MycoBank:  
      
    Vincent, R *et al* (2013). **MycoBank gearing up for new horizons.**
    IMA Fungus, 4(2), 371-9;
    [doi:10.5598/imafungus.2013.04.02.16](https://doi.org/10.5598/imafungus.2013.04.02.16)
    . Accessed from <https://www.mycobank.org> on 7th of May, 2026.

3.  Global Biodiversity Information Facility (GBIF), via Catalogue of
    Life (COL):  
      
    Banki, O. *et al.* (2026). Catalogue of Life (2026-04-18 XR).
    Catalogue of Life Foundation, Amsterdam, Netherlands.
    [doi:10.48580/dgxjw](https://doi.org/10.48580/dgxjw) . Accessed from
    <https://www.gbif.org> on 7th of May, 2026.

Furthermore, these sources were used for additional details:

- BacDive:  
    
  Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base
  for standardized bacterial and archaeal data.** Nucleic Acids Res.,
  50(D1):D741-D74;
  [doi:10.1093/nar/gkab961](https://doi.org/10.1093/nar/gkab961) .
  Accessed from <https://bacdive.dsmz.de> on 7th of May, 2026.

- Systematized Nomenclature of Medicine - Clinical Terms (SNOMED-CT):  
    
  Public Health Information Network Vocabulary Access and Distribution
  System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
  Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).
  Accessed from <https://www.cdc.gov/phin/php/phinvads/> on 16th of
  July, 2024.

- Grimont *et al.* (2007). Antigenic Formulae of the Salmonella
  Serovars, 9th Edition. WHO Collaborating Centre for Reference and
  Research on *Salmonella* (WHOCC-SALM).

- Bartlett *et al.* (2022). **A comprehensive list of bacterial
  pathogens infecting humans** *Microbiology* 168:001269;
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## See also

[`as.mo()`](https://amr-for-r.org/reference/as.mo.md),
[`mo_property()`](https://amr-for-r.org/reference/mo_property.md),
[microorganisms.groups](https://amr-for-r.org/reference/microorganisms.groups.md),
[microorganisms.codes](https://amr-for-r.org/reference/microorganisms.codes.md),
[intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md)

## Examples

``` r
microorganisms
#> # A tibble: 96,982 × 28
#>    mo               fullname     status domain kingdom phylum class order family
#>    <mo>             <chr>        <chr>  <chr>  <chr>   <chr>  <chr> <chr> <chr> 
#>  1 B_GRAMN          (unknown Gr… unkno… Bacte… (unkno… (unkn… (unk… (unk… (unkn…
#>  2 B_GRAMP          (unknown Gr… unkno… Bacte… (unkno… (unkn… (unk… (unk… (unkn…
#>  3 B_ANAER-NEG      (unknown an… unkno… Bacte… (unkno… (unkn… (unk… (unk… (unkn…
#>  4 B_ANAER-POS      (unknown an… unkno… Bacte… (unkno… (unkn… (unk… (unk… (unkn…
#>  5 B_ANAER          (unknown an… unkno… Bacte… (unkno… (unkn… (unk… (unk… (unkn…
#>  6 F_FUNGUS         (unknown fu… unkno… Fungi  Fungi   (unkn… (unk… (unk… (unkn…
#>  7   UNKNOWN        (unknown na… unkno… (unkn… (unkno… (unkn… (unk… (unk… (unkn…
#>  8 P_PROTOZOAN      (unknown pr… unkno… Proto… Protoz… (unkn… (unk… (unk… (unkn…
#>  9 F_YEAST          (unknown ye… unkno… Fungi  Fungi   (unkn… (unk… (unk… (unkn…
#> 10 B_[FAM]_ABDTBCTR Abditibacte… accep… Bacte… Pseudo… Abdit… Abdi… Abdi… Abdit…
#> # ℹ 96,972 more rows
#> # ℹ 19 more variables: genus <chr>, species <chr>, subspecies <chr>,
#> #   rank <chr>, ref <chr>, oxygen_tolerance <chr>, morphology <chr>,
#> #   source <chr>, lpsn <chr>, lpsn_parent <chr>, lpsn_renamed_to <chr>,
#> #   mycobank <chr>, mycobank_parent <chr>, mycobank_renamed_to <chr>,
#> #   gbif <chr>, gbif_parent <chr>, gbif_renamed_to <chr>, prevalence <dbl>,
#> #   snomed <list>
```

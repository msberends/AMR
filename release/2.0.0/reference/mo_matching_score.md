# Calculate the Matching Score for Microorganisms

This algorithm is used by
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md) and
all the
[`mo_*`](https://amr-for-r.org/release/2.0.0/reference/mo_property.md)
functions to determine the most probable match of taxonomic records
based on user input.

## Usage

``` r
mo_matching_score(x, n)
```

## Arguments

- x:

  Any user input value(s)

- n:

  A full taxonomic name, that exists in
  [`microorganisms$fullname`](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)

## Note

This algorithm was originally described in: Berends MS *et al.* (2022).
**AMR: An R Package for Working with Antimicrobial Resistance Data**.
*Journal of Statistical Software*, 104(3), 1-31;
[doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03) .

Later, the work of Bartlett A *et al.* about bacterial pathogens
infecting humans (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ) was
incorporated.

## Matching Score for Microorganisms

With ambiguous user input in
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md) and
all the
[`mo_*`](https://amr-for-r.org/release/2.0.0/reference/mo_property.md)
functions, the returned results are chosen based on their matching score
using `mo_matching_score()`. This matching score \\m\\, is calculated
as:

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
  = 2, Protozoa = 3, Archaea = 4, others = 5.

The grouping into human pathogenic prevalence \\p\\ is based on recent
work from Bartlett *et al.* (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ) who
extensively studied medical-scientific literature to categorise all
bacterial species into these groups:

- **Established**, if a taxonomic species has infected at least three
  persons in three or more references. These records have
  `prevalence = 1.0` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set;

- **Putative**, if a taxonomic species has fewer than three known cases.
  These records have `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set.

Furthermore,

- Any genus present in the **established** list also has
  `prevalence = 1.0` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set;

- Any other genus present in the **putative** list has
  `prevalence = 1.25` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set;

- Any other species or subspecies of which the genus is present in the
  two aforementioned groups, has `prevalence = 1.5` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set;

- Any *non-bacterial* genus, species or subspecies of which the genus is
  present in the following list, has `prevalence = 1.5` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
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
  *Phlebotomus*, *Phoma*, *Pichia*, *Piedraia*, *Pithomyces*,
  *Pityrosporum*, *Pneumocystis*, *Pseudallescheria*, *Pseudoterranova*,
  *Pulex*, *Rhizomucor*, *Rhizopus*, *Rhodotorula*, *Saccharomyces*,
  *Sarcoptes*, *Scolecobasidium*, *Scopulariopsis*, *Scytalidium*,
  *Spirometra*, *Sporobolomyces*, *Stachybotrys*, *Strongyloides*,
  *Syngamus*, *Taenia*, *Toxocara*, *Trichinella*, *Trichobilharzia*,
  *Trichoderma*, *Trichomonas*, *Trichophyton*, *Trichosporon*,
  *Trichostrongylus*, *Trichuris*, *Tritirachium*, *Trombicula*,
  *Trypanosoma*, *Tunga*, or *Wuchereria*;

- All other records have `prevalence = 2.0` in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set.

When calculating the matching score, all characters in \\x\\ and \\n\\
are ignored that are other than A-Z, a-z, 0-9, spaces and parentheses.

All matches are sorted descending on their matching score and for all
user input values, the top match will be returned. This will lead to the
effect that e.g., `"E. coli"` will return the microbial ID of
*Escherichia coli* (\\m = 0.688\\, a highly prevalent microorganism
found in humans) and not *Entamoeba coli* (\\m = 0.159\\, a less
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

## Author

Dr. Matthijs Berends, 2018

## Examples

``` r
mo_reset_session()

as.mo("E. coli")
mo_uncertainties()

mo_matching_score(
  x = "E. coli",
  n = c("Escherichia coli", "Entamoeba coli")
)
```

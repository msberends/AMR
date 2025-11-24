# Calculate the Matching Score for Microorganisms

This algorithm is used by
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all the
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions to
determine the most probable match of taxonomic records based on user
input.

## Usage

``` r
mo_matching_score(x, n)
```

## Arguments

- x:

  Any user input value(s).

- n:

  A full taxonomic name, that exists in
  [`microorganisms$fullname`](https://amr-for-r.org/reference/microorganisms.md).

## Note

This algorithm was originally developed in 2018 and subsequently
described in: Berends MS *et al.* (2022). **AMR: An R Package for
Working with Antimicrobial Resistance Data**. *Journal of Statistical
Software*, 104(3), 1-31;
[doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03) .

Later, the work of Bartlett A *et al.* about bacterial pathogens
infecting humans (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ) was
incorporated, and optimalisations to the algorithm were made.

## Matching Score for Microorganisms

With ambiguous user input in
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and all the
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions, the
returned results are chosen based on their matching score using
`mo_matching_score()`. This matching score \\m\\, is calculated as:

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

## Examples

``` r
mo_reset_session()
#> ℹ Reset 17 previously matched input values.

as.mo("E. coli")
#> Class 'mo'
#> [1] B_ESCHR_COLI
mo_uncertainties()
#> Matching scores are based on the resemblance between the input and the full
#> taxonomic name, and the pathogenicity in humans. See `?mo_matching_score`.
#> Colour keys:  0.000-0.549  0.550-0.649  0.650-0.749  0.750-1.000 
#> 
#> --------------------------------------------------------------------------------
#> "E. coli" -> Escherichia coli (B_ESCHR_COLI, 0.688)
#> Also matched: Enterococcus crotali (0.650), Escherichia coli coli
#>               (0.643), Escherichia coli expressing (0.611), Enterobacter cowanii
#>               (0.600), Enterococcus columbae (0.595), Enterococcus camelliae (0.591),
#>               Enterococcus casseliflavus (0.577), Enterobacter cloacae cloacae
#>               (0.571), Enterobacter cloacae complex (0.571), and Enterobacter cloacae
#>               dissolvens (0.565)
#> 
#> Only the first 10 other matches of each record are shown. Run
#> `print(mo_uncertainties(), n = ...)` to view more entries, or save
#> `mo_uncertainties()` to an object.

mo_matching_score(
  x = "E. coli",
  n = c("Escherichia coli", "Entamoeba coli")
)
#> [1] 0.6875000 0.3809524
```

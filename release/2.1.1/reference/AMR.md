# The `AMR` Package

Welcome to the `AMR` package.

The `AMR` package is a [free and
open-source](https://msberends.github.io/AMR/#copyright) R package with
[zero dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to
simplify the analysis and prediction of Antimicrobial Resistance (AMR)
and to work with microbial and antimicrobial data and properties, by
using evidence-based methods. **Our aim is to provide a standard** for
clean and reproducible AMR data analysis, that can therefore empower
epidemiological analyses to continuously enable surveillance and
treatment evaluation in any setting. [Many different
researchers](https://msberends.github.io/AMR/authors.html) from around
the globe are continually helping us to make this a successful and
durable project!

This work was published in the Journal of Statistical Software (Volume
104(3);
[doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03) ) and
formed the basis of two PhD theses
([doi:10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131)
and
[doi:10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)
).

After installing this package, R knows [**~52 000
microorganisms**](https://msberends.github.io/AMR/reference/microorganisms.html)
(updated December 2022) and all [**~600 antibiotic, antimycotic and
antiviral
drugs**](https://msberends.github.io/AMR/reference/antibiotics.html) by
name and code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and
SNOMED CT), and knows all about valid SIR and MIC values. The integral
clinical breakpoint guidelines from CLSI and EUCAST are included, even
with epidemiological cut-off (ECOFF) values. It supports and can read
any data format, including WHONET data. This package works on Windows,
macOS and Linux with all versions of R since R-3.0 (April 2013). **It
was designed to work in any setting, including those with very limited
resources**. It was created for both routine data analysis and academic
research at the Faculty of Medical Sciences of the public [University of
Groningen](https://www.rug.nl), in collaboration with non-profit
organisations [Certe Medical Diagnostics and Advice
Foundation](https://www.certe.nl) and [University Medical Center
Groningen](https://www.umcg.nl).

The `AMR` package is available in English, Chinese, Czech, Danish,
Dutch, Finnish, French, German, Greek, Italian, Japanese, Norwegian,
Polish, Portuguese, Romanian, Russian, Spanish, Swedish, Turkish, and
Ukrainian. Antimicrobial drug (group) names and colloquial microorganism
names are provided in these languages.

## Source

To cite AMR in publications use:

Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C
(2022). "AMR: An R Package for Working with Antimicrobial Resistance
Data." *Journal of Statistical Software*, *104*(3), 1-31.
[doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)

A BibTeX entry for LaTeX users is:


    @Article{,
      title = {{AMR}: An {R} Package for Working with Antimicrobial Resistance Data},
      author = {Matthijs S. Berends and Christian F. Luz and Alexander W. Friedrich and Bhanu N. M. Sinha and Casper J. Albers and Corinna Glasner},
      journal = {Journal of Statistical Software},
      year = {2022},
      volume = {104},
      number = {3},
      pages = {1--31},
      doi = {10.18637/jss.v104.i03},
    }

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

Useful links:

- <https://msberends.github.io/AMR/>

- <https://github.com/msberends/AMR>

- Report bugs at <https://github.com/msberends/AMR/issues>

## Author

**Maintainer**: Matthijs S. Berends <m.s.berends@umcg.nl>
([ORCID](https://orcid.org/0000-0001-7620-1800))

Authors:

- Christian F. Luz ([ORCID](https://orcid.org/0000-0001-5809-5995))
  \[contributor\]

- Dennis Souverein ([ORCID](https://orcid.org/0000-0003-0455-0336))
  \[contributor\]

- Erwin E. A. Hassing \[contributor\]

Other contributors:

- Casper J. Albers ([ORCID](https://orcid.org/0000-0002-9213-6743))
  \[thesis advisor\]

- Peter Dutey-Magni ([ORCID](https://orcid.org/0000-0002-8942-9836))
  \[contributor\]

- Judith M. Fonville \[contributor\]

- Alex W. Friedrich ([ORCID](https://orcid.org/0000-0003-4881-038X))
  \[thesis advisor\]

- Corinna Glasner ([ORCID](https://orcid.org/0000-0003-1241-1328))
  \[thesis advisor\]

- Eric H. L. C. M. Hazenberg \[contributor\]

- Gwen Knight ([ORCID](https://orcid.org/0000-0002-7263-9896))
  \[contributor\]

- Annick Lenglet ([ORCID](https://orcid.org/0000-0003-2013-8405))
  \[contributor\]

- Bart C. Meijer \[contributor\]

- Dmytro Mykhailenko \[contributor\]

- Anton Mymrikov \[contributor\]

- Andrew P. Norgan ([ORCID](https://orcid.org/0000-0002-2955-2066))
  \[contributor\]

- Sofia Ny ([ORCID](https://orcid.org/0000-0002-2017-1363))
  \[contributor\]

- Jonas Salm \[contributor\]

- Rogier P. Schade \[contributor\]

- Bhanu N. M. Sinha ([ORCID](https://orcid.org/0000-0003-1634-0010))
  \[thesis advisor\]

- Anthony Underwood ([ORCID](https://orcid.org/0000-0002-8547-4277))
  \[contributor\]

- Anita Williams ([ORCID](https://orcid.org/0000-0002-5295-8451))
  \[contributor\]

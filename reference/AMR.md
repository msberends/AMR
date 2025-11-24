# The `AMR` Package

Welcome to the `AMR` package.

The `AMR` package is a peer-reviewed, [free and
open-source](https://amr-for-r.org/#copyright) R package with [zero
dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to simplify
the analysis and prediction of Antimicrobial Resistance (AMR) and to
work with microbial and antimicrobial data and properties, by using
evidence-based methods. **Our aim is to provide a standard** for clean
and reproducible AMR data analysis, that can therefore empower
epidemiological analyses to continuously enable surveillance and
treatment evaluation in any setting. We are a team of [many different
researchers](https://amr-for-r.org/authors.html) from around the globe
to make this a successful and durable project!

This work was published in the Journal of Statistical Software (Volume
104(3);
[doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03) ) and
formed the basis of two PhD theses
([doi:10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131)
and
[doi:10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)
).

After installing this package, R knows [**~79 000 distinct microbial
species**](https://amr-for-r.org/reference/microorganisms.html) (updated
June 2024) and all [**~620 antimicrobial and antiviral
drugs**](https://amr-for-r.org/reference/antimicrobials.html) by name
and code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and SNOMED
CT), and knows all about valid SIR and MIC values. The integral clinical
breakpoint guidelines from CLSI 2011-2025 and EUCAST 2011-2025 are
included, even with epidemiological cut-off (ECOFF) values. It supports
and can read any data format, including WHONET data. This package works
on Windows, macOS and Linux with all versions of R since R-3.0 (April
2013). **It was designed to work in any setting, including those with
very limited resources**. It was created for both routine data analysis
and academic research at the Faculty of Medical Sciences of the
[University of Groningen](https://www.rug.nl) and the [University
Medical Center Groningen](https://www.umcg.nl).

The `AMR` package is available in English, Arabic, Bengali, Chinese,
Czech, Danish, Dutch, Finnish, French, German, Greek, Hindi, Indonesian,
Italian, Japanese, Korean, Norwegian, Polish, Portuguese, Romanian,
Russian, Spanish, Swahili, Swedish, Turkish, Ukrainian, Urdu, and
Vietnamese. Antimicrobial drug (group) names and colloquial
microorganism names are provided in these languages.

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

Useful links:

- <https://amr-for-r.org>

- <https://github.com/msberends/AMR>

- Report bugs at <https://github.com/msberends/AMR/issues>

## Author

**Maintainer**: Matthijs S. Berends <m.s.berends@umcg.nl>
([ORCID](https://orcid.org/0000-0001-7620-1800))

Authors:

- Dennis Souverein ([ORCID](https://orcid.org/0000-0003-0455-0336))
  \[contributor\]

- Erwin E. A. Hassing \[contributor\]

Other contributors:

- Aislinn Cook ([ORCID](https://orcid.org/0000-0002-9189-7815))
  \[contributor\]

- Andrew P. Norgan ([ORCID](https://orcid.org/0000-0002-2955-2066))
  \[contributor\]

- Anita Williams ([ORCID](https://orcid.org/0000-0002-5295-8451))
  \[contributor\]

- Annick Lenglet ([ORCID](https://orcid.org/0000-0003-2013-8405))
  \[contributor\]

- Anthony Underwood ([ORCID](https://orcid.org/0000-0002-8547-4277))
  \[contributor\]

- Anton Mymrikov \[contributor\]

- Bart C. Meijer \[contributor\]

- Christian F. Luz ([ORCID](https://orcid.org/0000-0001-5809-5995))
  \[contributor\]

- Dmytro Mykhailenko \[contributor\]

- Eric H. L. C. M. Hazenberg \[contributor\]

- Gwen Knight ([ORCID](https://orcid.org/0000-0002-7263-9896))
  \[contributor\]

- Jane Hawkey ([ORCID](https://orcid.org/0000-0001-9661-5293))
  \[contributor\]

- Jason Stull ([ORCID](https://orcid.org/0000-0002-9028-8153))
  \[contributor\]

- Javier Sanchez ([ORCID](https://orcid.org/0000-0003-2605-8094))
  \[contributor\]

- Jonas Salm \[contributor\]

- Judith M. Fonville \[contributor\]

- Kathryn Holt ([ORCID](https://orcid.org/0000-0003-3949-2471))
  \[contributor\]

- Larisse Bolton ([ORCID](https://orcid.org/0000-0001-7879-2173))
  \[contributor\]

- Matthew Saab ([ORCID](https://orcid.org/0009-0008-6626-7919))
  \[contributor\]

- Natacha Couto ([ORCID](https://orcid.org/0000-0002-9152-5464))
  \[contributor\]

- Peter Dutey-Magni ([ORCID](https://orcid.org/0000-0002-8942-9836))
  \[contributor\]

- Rogier P. Schade ([ORCID](https://orcid.org/0000-0002-9487-4467))
  \[contributor\]

- Sofia Ny ([ORCID](https://orcid.org/0000-0002-2017-1363))
  \[contributor\]

- Alex W. Friedrich ([ORCID](https://orcid.org/0000-0003-4881-038X))
  \[thesis advisor\]

- Bhanu N. M. Sinha ([ORCID](https://orcid.org/0000-0003-1634-0010))
  \[thesis advisor\]

- Casper J. Albers ([ORCID](https://orcid.org/0000-0002-9213-6743))
  \[thesis advisor\]

- Corinna Glasner ([ORCID](https://orcid.org/0000-0003-1241-1328))
  \[thesis advisor\]

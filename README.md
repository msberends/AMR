# `AMR` (for R)

<a href="https://msberends.github.io/AMR/"><img src="https://msberends.github.io/AMR/AMR_intro.png" align="center"></a>

----

This work was published in the Journal of Statistical Software (Volume 104(3); [DOI 10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)) and formed the basis of two PhD theses ([DOI 10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131) and [DOI 10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)).

`AMR` is a free, open-source and independent R package to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. Our aim is to provide a standard for clean and reproducible antimicrobial resistance data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting. It is currently being used in over 175 countries.
 
After installing this package, R knows ~52,000 distinct microbial species and all ~600 antibiotic, antimycotic, and antiviral drugs by name and code (including ATC, WHONET/EARS-Net, PubChem, LOINC and SNOMED CT), and knows all about valid SIR and MIC values. It supports any data format, including WHONET/EARS-Net data. Antimicrobial names and group names are available in English, Chinese, Danish, Dutch, French, German, Greek, Italian, Japanese, Polish, Portuguese, Russian, Spanish, Swedish, Turkish, and Ukrainian.

This package is fully independent of any other R package and works on Windows, macOS and Linux with all versions of R since R-3.0.0 (April 2013). It was designed to work in any setting, including those with very limited resources. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the University of Groningen, in collaboration with non-profit organisations Certe Medical Diagnostics and Advice Foundation and University Medical Center Groningen. This R package is actively maintained and free software; you can freely use and distribute it for both personal and commercial (but not patent) purposes under the terms of the GNU General Public License version 2.0 (GPL-2), as published by the Free Software Foundation.

This is the development source of the `AMR` package for R. Not a developer? Then please visit our website [https://msberends.github.io/AMR/](https://msberends.github.io/AMR/) to read more about this package.

*NOTE: this source code is on GitHub (https://github.com/msberends/AMR), but also automatically mirrored to our university's Gitea server (https://git.web.rug.nl/P281424/AMR) and to GitLab (https://gitlab.com/msberends/AMR).*

### How to get this package
Please see [our website](https://msberends.github.io/AMR/#get-this-package).

You can install or update the `AMR` package from CRAN using:

```r
install.packages("AMR")
```

It will be downloaded and installed automatically. For RStudio, click on the menu *Tools* > *Install Packages...* and then type in "AMR" and press <kbd>Install</kbd>.

### Copyright

This R package is licensed under the [GNU General Public License (GPL) v2.0](https://github.com/msberends/AMR/blob/main/LICENSE). In a nutshell, this means that this package:

- May be used for commercial purposes

- May be used for private purposes

- May **not** be used for patent purposes

- May be modified, although:

  - Modifications **must** be released under the same license when distributing the package
  - Changes made to the code **must** be documented

- May be distributed, although:

  - Source code **must** be made available when the package is distributed
  - A copy of the license and copyright notice **must** be included with the package.

- Comes with a LIMITATION of liability

- Comes with NO warranty

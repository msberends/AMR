# `AMR` (for R) <img src="man/figures/logo.png" align="right" height="120px" />

*NOTE: the original source code is on GitLab (https://gitlab.com/msberends/AMR), so you can report a bug at https://gitlab.com/msberends/AMR/issues. There is a mirror repository on GitHub (https://github.com/msberends/AMR). As the mirror process is automated by GitLab, both repositories always contain the latest changes.*

----

## Development source

This is the **development source** of `AMR`, a free and open-source [R package](https://www.r-project.org) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with antibiotic properties by using evidence-based methods.

**Not a developer? Then our website https://msberends.gitlab.io/AMR is probably a better place to read about this package.** It contains documentation about all of the included functions and also a comprehensive tutorial about how to conduct AMR analysis.

## Authors
Matthijs S. Berends <a href="https://orcid.org/0000-0001-7620-1800"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,2,a</sup>,
Christian F. Luz <a href="https://orcid.org/0000-0001-5809-5995"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,a</sup>,
Erwin E.A. Hassing<sup>2</sup>,
Corinna Glasner <a href="https://orcid.org/0000-0003-1241-1328"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>,
Alex W. Friedrich <a href="https://orcid.org/0000-0003-4881-038X"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>,
Bhanu N.M. Sinha <a href="https://orcid.org/0000-0003-1634-0010"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>
  
<sup>1</sup> Department of Medical Microbiology, University of Groningen, University Medical Center Groningen, Groningen, the Netherlands - [rug.nl](http://www.rug.nl) [umcg.nl](http://www.umcg.nl)<br>
<sup>2</sup> Certe Medical Diagnostics & Advice, Groningen, the Netherlands - [certe.nl](http://www.certe.nl)<br>
<sup>a</sup> R package author and thesis dissertant<br>
<sup>b</sup> Thesis advisor

<a href="https://www.rug.nl"><img src="man/figures/logo_rug.png" height="60px"></a>
<a href="https://www.umcg.nl"><img src="man/figures/logo_umcg.png" height="60px"></a>
<a href="https://www.certe.nl"><img src="man/figures/logo_certe.png" height="60px"></a>
<a href="http://www.eurhealth-1health.eu"><img src="man/figures/logo_eh1h.png" height="60px"></a>
<a href="http://www.eurhealth-1health.eu"><img src="man/figures/logo_interreg.png" height="60px"></a>

## How to get this package
All stable versions of this package [are published on CRAN](https://CRAN.R-project.org/package=AMR), the official R network with a peer-reviewed submission process.

### Install from CRAN
[![CRAN_Badge](https://www.r-pkg.org/badges/version/AMR)](https://CRAN.R-project.org/package=AMR) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/AMR)](https://CRAN.R-project.org/package=AMR)

(Note: Downloads measured only by [cran.rstudio.com](https://cran.rstudio.com/package=AMR), this excludes e.g. the official [cran.r-project.org](https://cran.r-project.org/package=AMR))

- <img src="http://www.rstudio.com/favicon.ico" alt="RStudio favicon" height="20px"> Install using [RStudio](http://www.rstudio.com) (recommended):
  - Click on `Tools` and then `Install Packages...`
  - Type in `AMR` and press <kbd>Install</kbd>

- <img src="https://cran.r-project.org/favicon.ico" alt="R favicon" height="20px"> Install in R directly:
  - `install.packages("AMR")`

### Install from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1305355.svg)](https://doi.org/10.5281/zenodo.1305355)

This package was also published on Zenodo (stable releases only): https://doi.org/10.5281/zenodo.1305355

### Install from GitLab

This is the latest **development version**. Although it may contain bugfixes and even new functions compared to the latest released version on CRAN, it is also subject to change and may be unstable or behave unexpectedly. Always consider this a beta version. All below 'badges' should be green:

Development Test | Result | Reference
--- | :---: | ---
All functions checked on Linux | [![pipeline status](https://gitlab.com/msberends/AMR/badges/master/pipeline.svg)](https://gitlab.com/msberends/AMR/commits/master) | GitLab CI [[ref 1]](https://gitlab.com/msberends/AMR/pipelines) 
All functions checked on Windows | [![AppVeyor_Build](https://ci.appveyor.com/api/projects/status/gitlab/msberends/AMR?branch=master&svg=true)](https://ci.appveyor.com/project/msberends/amr-svxon) | Appveyor Systems Inc. [[ref 2]](https://ci.appveyor.com/project/msberends/amr-svxon)
Percentage of syntax lines checked | [![Code_Coverage](https://codecov.io/gl/msberends/AMR/branch/master/graph/badge.svg)](https://codecov.io/gl/msberends/AMR) | Codecov LLC [[ref 3]](https://codecov.io/gl/msberends/AMR)

If so, try it with:
```r
install.packages("devtools") 
devtools::install_gitlab("msberends/AMR")
```

## Copyright

This R package is licensed under the [GNU General Public License (GPL) v2.0](https://gitlab.com/msberends/AMR/blob/master/LICENSE). In a nutshell, this means that this package:

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

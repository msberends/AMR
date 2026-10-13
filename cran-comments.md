This version is a minor release (v3.1.0) following the bugfix release v3.0.1 in September 2025. It adds the EUCAST 2026 and CLSI 2026 clinical breakpoints (EUCAST 2026 is now the default), a taxonomic update of all microorganisms, unified EUCAST and CLSI interpretive rules, support for `tidymodels` pipelines, and many improvements and bug fixes, as listed in `NEWS.md`.

This package has two reverse dependencies, `MIC` and `MIMER`, both checked against this version. `MIMER` passes. `MIC` is affected by a change in the output of `as.sir()` for interpretation based on epidemiological cut-off values (ECOFFs), which makes some of its unit tests fail; its maintainer has been informed in advance through an issue in its GitHub repository.

As with all previous >20 releases, some CHECKs on `oldrel` may return a `NOTE` for narrowly exceeding the installation size limit. This has been reduced to a minimum in prior coordination with CRAN maintainers and currently returns only an `INFO` on `release` and `devel`.

We treat this as a high-impact package: it was published in the *Journal of Statistical Software* (2022), is listed in the CRAN Task View "Epidemiology", and (based on cranlogs download statistics) is used globally. If there is anything to address, we would appreciate being informed before archiving the current version. We conduct extensive automated unit testing and have no indication of unresolved issues.

Thank you for your continued maintenance of CRAN, it plays a central role in the success and growth of the R ecosystem.

This version is a minor release (v3.1.0) following the bugfix release v3.0.1 in September 2025.The many improvements and bug fixes are in `NEWS.md`.

This package has two reverse dependencies, `MIC` and `MIMER`, both checked against this version. `MIMER` passes. `MIC` is affected by a change in the output of `AMR::as.sir()` for interpretation based on epidemiological cut-off values (ECOFFs), which makes some of its unit tests fail; its maintainer has been informed [in advance through this issue](https://github.com/agerada/MIC/issues/15) in their GitHub repository.

As with all previous >20 releases, some CHECKs may return an `INFO` for narrowly exceeding the installation size limit. The size has been reduced to a minimum in prior coordination with CRAN maintainers.

We treat this as a high-impact package: it was published in the *Journal of Statistical Software* (2022), is listed in the CRAN Task View "Epidemiology", and (based on cranlogs download statistics) is used globally. If there is anything to address, we would appreciate being informed before archiving the current version. We conduct extensive automated unit testing and have no indication of unresolved issues.

Thank you for your continued maintenance of CRAN, it plays a central role in the success and growth of the R ecosystem.

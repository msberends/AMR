* Edited the unit tests, so they will run under 10 minutes on CRAN (using testthat::skip_on_cran() on some tests).

* Since version 0.3.0 (2018-08-14), CHECK returns a NOTE for having a data directory over 3 MB. This is needed to offer users reference data for the complete taxonomy of microorganisms - one of the most important features of this package.

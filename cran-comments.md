# Version 0.8.0

* A NOTE for having a data directory over 3 MB. This is needed to offer users reference data for the complete taxonomy of microorganisms - one of the most important features of this pacakge. Has been this way since version 0.3.0.

* This package writes lines to `[library path]/AMR/mo_history/mo_history.csv` when using the `as.mo()` function, in the exact same way (and borrowed from) the `extrafont` package on CRAN (version 0.17) writes to the user library path. Users are notified about this with a `message()`, and staged install on R >= 3.6.0 still works. The CSV file is never newly created or deleted by this package, it only changes this file to improve speed and reliability of the `as.mo()` function. See the source code of functions `set_mo_history()` and `clear_mo_history()` in file `R/mo_history.R`.

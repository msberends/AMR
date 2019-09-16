# Version 0.8.0

* A NOTE for having a data directory over 3 MB. This is needed to offer users reference data for the complete taxonomy of microorganisms - one of the most important features of this pacakge. Has been this way since version 0.3.0.
* This package writes lines to `[user library]/AMR/inst/mo_history/mo_history.csv` when using the `as.mo()` function. Users are notified about this. The CSV file is never newly created or deleted by this package, it only changes this file to improve speed and reliability of the `as.mo()` function. Staged install still works. The source code was taken from the `extrafont` package on CRAN (version 0.17), that writes to the package folder in the user library exactly the same way. See the source code of `set_mo_history()` and `clear_mo_history()`.

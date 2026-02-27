# CLAUDE.md — AMR R Package

This file provides context for Claude Code when working in this
repository.

## Project Overview

**AMR** is a zero-dependency R package for antimicrobial resistance
(AMR) data analysis using a One Health approach. It is peer-reviewed,
used in 175+ countries, and supports 28 languages.

Key capabilities: - SIR (Susceptible/Intermediate/Resistant)
classification using EUCAST 2011–2025 and CLSI 2011–2025 breakpoints -
Antibiogram generation: traditional, combined, syndromic, and WISCA -
Microorganism taxonomy database (~79,000 species) - Antimicrobial drug
database (~620 drugs) - Multi-drug resistant organism (MDRO)
classification - First-isolate identification - Minimum Inhibitory
Concentration (MIC) and disk diffusion handling - Multilingual output
(28 languages)

## Common Commands

All commands run inside an R session:

``` r
# Rebuild documentation (roxygen2 → .Rd files + NAMESPACE)
devtools::document()

# Run all tests
devtools::test()

# Full package check (CRAN-level: docs + tests + checks)
devtools::check()

# Build pkgdown website locally
pkgdown::build_site()

# Code coverage report
covr::package_coverage()
```

From the shell:

``` bash
# CRAN check from parent directory
R CMD check AMR
```

## Repository Structure

    R/              # All R source files (62 files, ~28,000 lines)
    man/            # Auto-generated .Rd documentation (do not edit manually)
    tests/testthat/ # testthat test files (test-*.R) and helper-functions.R
    data/           # Pre-compiled .rda datasets
    data-raw/       # Scripts used to generate data/ files
    vignettes/      # Rmd vignette articles
    inst/           # Installed files (translations, etc.)
    _pkgdown.yml    # pkgdown website configuration

## R Source File Conventions

**Naming conventions in `R/`:**

| Prefix/Name       | Purpose                                                |
|-------------------|--------------------------------------------------------|
| `aa_*.R`          | Loaded first (helpers, globals, options, package docs) |
| `zz_deprecated.R` | Deprecated function wrappers                           |
| `zzz.R`           | `.onLoad` / `.onAttach` initialization                 |

**Key source files:**

- `aa_helper_functions.R` / `aa_helper_pm_functions.R` — internal
  utility functions (large; ~63 KB and ~37 KB)
- `aa_globals.R` — global constants and breakpoint lookup structures
- `aa_options.R` — `amr_options()` / `get_AMR_option()` system
- `mo.R` / `mo_property.R` — microorganism lookup and properties
- `ab.R` / `ab_property.R` — antimicrobial drug functions
- `av.R` / `av_property.R` — antiviral drug functions
- `sir.R` / `sir_calc.R` / `sir_df.R` — SIR classification engine
- `mic.R` / `disk.R` — MIC and disk diffusion classes
- `antibiogram.R` — antibiogram generation (traditional, combined,
  syndromic, WISCA)
- `first_isolate.R` — first-isolate identification algorithms
- `mdro.R` — MDRO classification (EUCAST, CLSI, CDC, custom guidelines)
- `amr_selectors.R` — tidyselect helpers for selecting AMR columns
- `interpretive_rules.R` / `custom_eucast_rules.R` — clinical
  interpretation rules
- `translate.R` — 28-language translation system
- `ggplot_sir.R` / `ggplot_pca.R` / `plotting.R` — visualisation
  functions

## Custom S3 Classes

The package defines five S3 classes with full print/format/plot/vctrs
support:

| Class    | Created by                                                | Represents                       |
|----------|-----------------------------------------------------------|----------------------------------|
| `<mo>`   | [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)     | Microorganism code               |
| `<ab>`   | [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)     | Antimicrobial drug code          |
| `<av>`   | [`as.av()`](https://amr-for-r.org/reference/as.av.md)     | Antiviral drug code              |
| `<sir>`  | [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)   | SIR value (S/I/R/SDD)            |
| `<mic>`  | [`as.mic()`](https://amr-for-r.org/reference/as.mic.md)   | Minimum inhibitory concentration |
| `<disk>` | [`as.disk()`](https://amr-for-r.org/reference/as.disk.md) | Disk diffusion diameter          |

## Data Files

Pre-compiled in `data/` (do not edit directly; regenerate via
`data-raw/` scripts):

| File                       | Contents                                      |
|----------------------------|-----------------------------------------------|
| `microorganisms.rda`       | ~79,000 microbial species with full taxonomy  |
| `antimicrobials.rda`       | ~620 antimicrobial drugs with ATC codes       |
| `antivirals.rda`           | Antiviral drugs                               |
| `clinical_breakpoints.rda` | EUCAST + CLSI breakpoints (2011–2025)         |
| `intrinsic_resistant.rda`  | Intrinsic resistance patterns                 |
| `example_isolates.rda`     | Example AMR dataset for documentation/testing |
| `WHONET.rda`               | Example WHONET-format dataset                 |

## Zero-Dependency Design

The package has **no `Imports`** in `DESCRIPTION`. All optional
integrations (ggplot2, dplyr, data.table, tidymodels, cli, crayon, etc.)
are listed in `Suggests` and guarded with:

``` r
if (requireNamespace("pkg", quietly = TRUE)) { ... }
```

Never add packages to `Imports`. If new functionality requires an
external package, add it to `Suggests` and guard usage appropriately.

## Testing

- **Framework:** `testthat` (R ≥ 3.1); legacy `tinytest` used for R
  3.0–3.6 CI
- **Test files:** `tests/testthat/test-*.R`
- **Helpers:** `tests/testthat/helper-functions.R`
- **CI matrix:** GitHub Actions across Windows / macOS / Linux × R devel
  / release / oldrel-1 through oldrel-4
- **Coverage:** `covr` (some files excluded: `atc_online.R`,
  `mo_source.R`, `translate.R`, `resistance_predict.R`,
  `zz_deprecated.R`, helper files, `zzz.R`)

## Documentation

- All exported functions use **roxygen2** blocks (`RoxygenNote: 7.3.3`,
  markdown enabled)
- Run `devtools::document()` after any change to roxygen comments
- Never edit files in `man/` directly — they are auto-generated
- Vignettes live in `vignettes/` as `.Rmd` files
- The pkgdown website is configured in `_pkgdown.yml`

## Versioning

Version format: `major.minor.patch.dev` (e.g., `3.0.1.9021`)

- Development versions use a `.9xxx` suffix
- Stable CRAN releases drop the dev suffix (e.g., `3.0.1`)
- `NEWS.md` uses sections **New**, **Fixes**, **Updates** with GitHub
  issue references (`#NNN`)

### Version bump required for every PR

Before opening a pull request, always increment the four-digit dev
counter by 1 in **both** of these files:

1.  **`DESCRIPTION`** — the `Version:` field:

        Version: 3.0.1.9021  →  Version: 3.0.1.9022

2.  **`NEWS.md`** — the top-level heading:

        # AMR 3.0.1.9021  →  # AMR 3.0.1.9022

Read the current version from `DESCRIPTION`, add 1 to the last numeric
component, and write the new version to both files in the same commit as
the rest of the PR changes.

## Internal State

The package uses a private `AMR_env` environment (created in
`aa_globals.R`) for caching expensive lookups (e.g., microorganism
matching scores, breakpoint tables). This avoids re-computation within a
session.

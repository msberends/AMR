# MO code registry

MO codes (such as `B_ESCHR_COLI`) are stored by users in their own data, sometimes for years. A code must
therefore **always mean the same taxon**, in every future version of the AMR package. This folder contains
the registry that guarantees this.

## Files

| File | Contents |
|---|---|
| `mo_code_registry.csv` | Every MO code that was ever part of a **released** version of the AMR package since v2.0.0, with the taxon it denotes (name, rank, domain, genus, species, subspecies), the first and last release it was in, and older names of the same code (e.g. a corrected spelling). |
| `build_mo_code_registry.R` | Recreates `mo_code_registry.csv` from the git tags of all releases. Run it after every new release. |
| `mo_code_renames.csv` | The only allowed exceptions: registered codes that get another name, because it is the same taxon (e.g. a corrected spelling or a renamed species group). Every row needs `approved_by` and `approved_date`. |

## Rules

1. **Only releases count.** Development versions (the main branch) are never part of the registry. Codes that
   only existed in a development version carry no weight and may be reassigned.
2. **Codes from before v2.0.0 are not supported.** Decision by Matthijs S. Berends, 3 October 2026: MO codes
   from releases before v2.0.0 (12 March 2023) are not supported. Every code in all releases since v2.0.0 has
   always meant the same taxon (checked on 3 October 2026), so the registry starts without conflicts.
3. **A registered code is never given to another taxon.** The taxonomy build
   (`data-raw/_reproduction_scripts/reproduction_of_microorganisms.R`) reuses the code of every taxon in the
   registry (most recent release first), never gives a registered code to a new taxon, and stops if a
   registered code would get another name than in the registry, unless that is listed and approved in
   `mo_code_renames.csv`.
4. **A registered code that is no longer in the data still works.** `as.mo()` translates it to the current
   code of the same taxon (by its name, so also if the taxon moved to another domain), using an internal lookup
   table (`MO_RETIRED_CODES`) that is created from the registry in `data-raw/_pre_commit_checks.R`. If the
   taxon is no longer in the data at all, the result is `NA` with a warning. The package never guesses.
5. **Tests enforce all of this.** `tests/testthat/test-data-microorganisms.R` checks the current data against
   the registry and the renames, and checks that every registered code resolves to its own taxon.

## After a release

1. Make sure the git tag of the release exists (e.g. `v3.1.0`).
2. Run `source("data-raw/microorganisms_files/build_mo_code_registry.R")` from the root of the repository.
3. Review the changes to `mo_code_registry.csv` (`git diff`): only new codes may be added, and existing codes
   may only get a later `last_release`. Commit the result.

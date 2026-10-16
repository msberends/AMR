# Release checklist

The steps of a regular CRAN release of the `AMR` package, as followed
for v3.1.0 (October 2026). Replace `X.Y.Z` with the new version number.
Copy this list into an issue or a note per release, or tick the boxes in
a local copy.

## 1. Preparation (on the development version)

All intended issues and pull requests for this release are merged into
`main`

All GitHub Actions on `main` succeed, including `check-revdeps` (it runs
on every push to `main` and can also be started manually from the
Actions tab)

Reverse dependencies are checked against this version; if a package
breaks, its maintainer is informed in advance through an issue in its
repository (e.g., <https://github.com/agerada/MIC/issues/15> for v3.1.0)

`vignettes/welcome_to_AMR.Rmd` is reviewed, since it is the only
vignette that goes to CRAN (all other vignettes are in `.Rbuildignore`):
check the text, the links, and the outcome of its inline R code (number
of species and drugs, guideline years, languages, taxonomy date)

`index.Rmd` (website home page) and `README.Rmd` are up to date with the
new version

`inst/CITATION` is up to date

## 2. Version and documentation

`DESCRIPTION`: `Version: X.Y.Z` and `Date:` set to the submission date

`NEWS.md`: line 1 replaced by `# AMR X.Y.Z`, the line about the planned
release removed, and the intro paragraph written in past tense,
summarising the main changes (see *NEWS.md style* in `CLAUDE.md`)

`NEWS.md`: the *Breaking* section mentions every change that can alter
results or break code, with how to get the old behaviour

`cran-comments.md` rewritten for this release: type of release, main
changes, reverse dependency results (with links to the issues of
affected packages), and the note on the installation size

`devtools::document()` is run, so that `man/`, `NAMESPACE` and
`R/sysdata.rda` are current (the pre-commit hook does this on commits
with the hooks active)

## 3. Checks

`R CMD build` and `R CMD check --as-cran` on the tarball: no errors or
warnings, and only the known note on the installation size

`devtools::check_win_devel()` (win-builder, R-devel); wait for the
e-mail

Optionally R-hub (`rhub::rhub_check()`) for other platforms

`urlchecker::url_check()` reports no broken or redirected URLs

## 4. Submission

`devtools::submit_cran()`; this also updates `CRAN-SUBMISSION` (version,
date and SHA)

Confirm the submission through the link in the e-mail from CRAN

Wait for acceptance; answer any questions from the CRAN team and
resubmit if needed (keep the version number unless CRAN asks otherwise)

## 5. After acceptance by CRAN

The pre-commit hook in `.github/prehooks/` sets the version to
`X.Y.Z.9NNN` and rewrites line 1 of `NEWS.md` on every commit to `main`.
The release commit must therefore skip the hooks.

`git status`: only the intended files are changed, with no build
leftovers (e.g., `AMR_X.Y.Z.tar.gz`, `AMR.Rcheck/`)

Commit the release without hooks, with only the version as message:
`git add -u && git commit -n -m "vX.Y.Z"`

`git push origin main` (this also rebuilds the website and publishes the
Python package to PyPI)

Create the tag through a GitHub release, with the `NEWS.md` section of
this version without its `# AMR X.Y.Z` heading as release notes:

``` bash
awk 'NR==1 {next} /^# AMR /{exit} 1' NEWS.md > notes.md
gh release create vX.Y.Z --target main --title "AMR vX.Y.Z" --notes-file notes.md --latest
rm notes.md
```

`git fetch --tags` and check that `vX.Y.Z` points at the release commit

Run `source("data-raw/microorganisms_files/build_mo_code_registry.R")`
from the root of the repository (takes a few minutes; the tag must
exist), then review the diff of `mo_code_registry.csv`: only new codes
may be added, existing codes may only get a later `last_release`, and
renamed codes must have their old name in `previous_names` (see
`data-raw/microorganisms_files/README.md`)

`.github/ISSUE_TEMPLATE/1-bug-report.yml`: set the versions to
`Latest CRAN version (X.Y.Z)` and
`One of the latest GitHub versions (X.Y.Z.9xxx)`

`NEWS.md`: put a new line 1 `# AMR X.Y.Z.9001`, followed by an empty
line, above `# AMR X.Y.Z`; without this, the hook overwrites the heading
of the released section

Commit with the hooks active (the version becomes `X.Y.Z.9001`) and push
to `main`

Check that the GitHub Actions succeed, the website
(<https://amr-for-r.org>) shows the new version, and the package is
available on PyPI

Check the CRAN page (<https://cran.r-project.org/package=AMR>) and the
CRAN check results some days later

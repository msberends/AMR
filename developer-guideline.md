# Developer Guideline

Welcome to the Developer Guideline of the `AMR` R package. This guideline explains about repository workflows and updates of package elements.

## Copyright

This R package and of its components are licensed under the [GNU General Public License (GPL) v2.0](https://github.com/msberends/AMR/blob/main/LICENSE). In a nutshell, this means that this package:

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

## General Intended Git(Hub) Workflow

All updates to the reposo should be done using `git commit`, preferably with the following predefined pre-commit Git hook.

This repository provides automated semantic versioning and R documentation updates by using a pre-commit Git hook. When using `git commit`, a script will be run to increase the version number, update the date and R documentation. To set this up, run this command when working locally in the repository:

```bash
git config --local core.hooksPath ".github/prehooks"
```

Now, when using `git commit`:

```bash
git commit -am "test commit"
# Running pre-commit hook...
# >>  Updating R documentation...
# >>  done.
# >>  
# >>  Updating semantic versioning and date...
# >>  - latest tag is 'v1.8.1', with 26 previous commits
# >>  - testpkg pkg version set to 1.8.1.9027
# >>  - updated DESCRIPTION
# >>  - updated NEWS.md
# >>  
# [main 300b93e] test commit
#  3 files changed, 3 insertions(+), 4 deletions(-)
```

### Website Generation

### Commiting Changes

## Updating the AMR Package

### Update EUCAST/CLSI Guidelines

### Update the Microbial Taxonomy

### Update the Antimicrobial Agents

### Add or Update a Language for Translation

For the most ideal workflow, please fork this repository and make the changes in your own forked repository. Afterwards, please create a Pull Request. If you are unfamiliar with these terms, no problem at all! Then please send us the edited files by email or any way you prefer.

The repository file [`data-raw/translations.tsv`](https://github.com/msberends/AMR/blob/main/developer-guideline.md) contains all translations. This file will be read by all functions where a translated output can be desired, like all `mo_*` functions (such as `mo_name()`, `mo_gramstain()`, `mo_type()`, etc.) and ``ab_*` functions (such as `ab_name()`, `ab_group()`, etc.). 

1.  To **add** a translation, edit `data-raw/translations.tsv` (you can copy the contains to MS Excel for convenience and paste the contents back later), add a column where the new column name is a [ISO 639-1 language code](https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes) (such as `en` for English, `de` for German and `es` for Spanish) and put in the new column all translated text from the first column.\
    \
    To **update** a translation, open `data-raw/translations.tsv` and save it with the language updates.

2.  Set the current working directory to the AMR package root folder (either by opening the AMR package as RStudio project, or by setting the working directory with `setwd()`).

3.  Run `source("data-raw/_language_update.R)"` to update the internal package data. If you have the `roxygen2` package installed, this script automatically updates the package documentation as well.

Many thanks for your contribution!


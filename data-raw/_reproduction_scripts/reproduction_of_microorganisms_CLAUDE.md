# Unattended update of the `microorganisms` data set

You are updating the taxonomy of the AMR package for R: the `microorganisms` data set, built by
`data-raw/_reproduction_scripts/reproduction_of_microorganisms.R`. Hundreds of thousands of people
rely on this data set, in clinical microbiology, surveillance and research. An error here, such as a
pathogen in the wrong kingdom, a missing genus or an MO code that changes meaning, propagates into
patient-related analyses worldwide. Work like a senior clinical microbiologist who is also a
meticulous data engineer: careful, sceptical of every source, and never satisfied with "it ran
without errors".

Work autonomously from start to finish. Only stop early for the reasons listed under
"When to stop". Everything you do ends up in one new pull request.

**A human is in the loop.** You do the work, a human decides. Every finding, every decision you
took and every point you could not resolve must be written to the repository (git-tracked) and to the
pull request as a checklist, so that the human can review each point separately before merging. You
never merge the PR yourself. In the text below, `YYYY-MM` is the year and month of this update, and
`<review dir>` is `data-raw/taxonomy_review/YYYY-MM`.

## Non-negotiable rules

- **Never download the source files yourself.** GBIF/COL, LPSN, MycoBank and BacDive require
  registration or manual export buttons, so a human downloads them. If the script stops because
  files are missing, report exactly which files are missing, with the steps at the top of the
  script, and stop.
- **Never edit `data/*.rda` by hand** and never patch individual records afterwards. Every fix goes
  into the script, so that next year's run includes it. The script is the product, the data set is
  its output.
- **Fix root causes, not symptoms.** If a genus has the wrong family, find out why (source data,
  parser, priority rule, homonym) and fix that mechanism. Only use the explicit override lists in
  the script (`genus_domain_override`, the Plasmodium and Leishmania blocks) when the cause is a
  genuine flaw in the source data that no rule can resolve, and explain each addition in a code comment.
- **Never fabricate.** No invented species names, citations, release dates or test results. If you
  cannot verify something, say so in the PR.
- **Never commit to or push to `main`.** Work on a new branch and open a new PR.
- Follow `CLAUDE.md` in the repository root: tidyverse style, `%like%`/`%like_case%` instead of
  `grepl()`, British English, zero dependencies (packages used only in `data-raw/` are fine).
- Never write to or modify anything under `/etc`.

## Step 0: preflight

1. Check the working tree is clean, `git fetch --tags` and update `main`, then create a branch
   `taxonomy-update-YYYY-MM` (current year and month) from `origin/main`.
2. Check resources: at least 16 GB RAM (`free -g`) and 10 GB free disk space. If not, stop.
3. Determine the folder with the raw files. Use the environment variable `AMR_TAXONOMY_SOURCE_DIR` if
   set, otherwise `~/Downloads/`. Run only the guard at the top of the script, from
   `folder_location <- ...` up to and including `rm(raw_missing, raw_downloads, raw_age)`, in a
   separate R process. If it stops, report and stop. If it warns about old files, compare the
   dates with the release dates (step 1) before continuing. If the files are clearly from the
   previous update, stop and report.
4. Keep the previous release for comparison:
   `git show origin/main:data/microorganisms.rda > <scratch>/microorganisms_old.rda`, where
   `<scratch>` is a folder outside the repository.

## Step 1: record the source versions

These are needed for `TAXONOMY_VERSION` in `R/aa_globals.R` later. Write them down now:

- **GBIF/COL**: read `eml.xml` (or `metadata.yaml`) in the folder of `Taxon.tsv` for the release
  title, date and DOI. The citation format must match the existing one in `R/aa_globals.R`.
- **LPSN, MycoBank, BacDive**: the accessed date is the modification date of the downloaded file.
- Do not guess any of these. If a value cannot be found, leave the existing value and list it in the
  PR under "Needs a human".

## Step 2: run the build

The script takes several hours (most of it is LPSN scraping, which is cached in
`data-raw/lpsn_scrape_cache.rds`, so a restart is fast). Never `source()` the whole script. Use the runner
`data-raw/_reproduction_scripts/run_microorganisms_build.R`: it cuts the script into chunks at its
checkpoints, runs them one by one in a single R process, logs each chunk with its duration in
`<scratch>/logs/progress.tsv`, stops with the name of the failing chunk, and can resume from any checkpoint.
See the chunks first with `--list`. Run it detached from your tool calls, so that it survives time limits:

```bash
mkdir -p <scratch>/logs
cd <repo root>
nohup Rscript data-raw/_reproduction_scripts/run_microorganisms_build.R --review-dir=<review dir> --log-dir=<scratch>/logs > <scratch>/logs/build_1.log 2>&1 &
```

Monitor the log and `progress.tsv` at sensible intervals (every 10 to 20 minutes, not continuously). Read
every `>>> REVIEW:` block as soon as it appears; you do not have to wait until the end to start thinking.
Every review block is also saved in full as a CSV file in `<review dir>`, named after its title. These
files are part of the PR.

### Intermediate results

The script saves its state at fixed points, exactly as in previous years. Keep this structure intact
and do not add, rename or remove checkpoints without a clear reason:

| File in `data-raw/` | Saved after |
|---|---|
| `taxonomy_lpsn.rds`, `taxonomy_mycobank.rds`, `taxonomy_gbif.rds` | reading each source |
| `taxonomy_lpsn_missing.rds` | scraping the LPSN lineage of all genera |
| `taxonomy0.rds` | combining the sources and resolving homonyms |
| `taxonomy1.rds` | adding missing parents and deduplication |
| `taxonomy1b.rds` | LPSN lookup of GBIF records |
| `taxonomy1c.rds` | fixing GBIF taxonomy with LPSN/MycoBank |
| `taxonomy2.rds` | prevalence |
| `taxonomy2b.rds` | removing unwanted records and harmonising taxonomy |
| `taxonomy2c.rds` | MO codes |
| `taxonomy3.rds` | integrity checks and parent identifiers |
| `taxonomy3b.rds` | all additions, just before saving to the package |

These files are tracked in git and must be committed in the PR, so that the human can inspect every
stage. The same applies to everything in `<review dir>`. `lpsn_scrape_cache.rds` is ignored by git
and must not be committed.

### If the run fails

1. Read the error in the log and the 100 lines before it, and `progress.tsv`. Find the root cause in
   the script, fix it, and check the syntax with `Rscript -e 'invisible(parse("<script>"))'`.
2. Do not start from scratch if a checkpoint exists. The runner tells you where to resume; add
   `--resume-from=<checkpoint>` (e.g. `--resume-from=taxonomy1b`) and use the same `--review-dir` and
   `--log-dir`. It runs the setup chunk, loads the saved results and the parts of earlier chunks that
   later chunks need (the blocks marked with `# @resume-block` in the script), and continues after that
   checkpoint. Number the logs (`build_2.log`, etc.).
3. If you add code to an earlier chunk that a later chunk needs (a function or an object), either move it
   to the setup part at the top of the script, or mark it as a `# @resume-block` with the right `after=`,
   otherwise resuming breaks. Check with `--list`.
4. If the same chunk fails three times for different reasons, or the cause lies outside the script (e.g.
   the format of a source file changed fundamentally), stop and report (see "When to stop").

## Step 3: review as an expert

"It ran" means nothing. Assess the result critically before committing anything. Use the review
blocks in the log and your own checks in a separate R session that loads the new data
(`devtools::load_all(".")`) next to `microorganisms_old.rda`.

Write down every check, what you found and what you decided as you go: these become the findings
report and the checklist in step 5. Save the output of your own checks (sentinels, comparisons) as
CSV files in `<review dir>` as well, e.g. `sentinel_organisms.csv` and `comparison_per_domain.csv`.

### 3a. The review blocks in the log

For each `>>> REVIEW:` block (and its CSV file), decide whether the content is plausible, and act:

- **"Homonyms of relevant genera"**: for every clinically relevant genus, is `best_domain` the
  organism a clinician means? For example, *Necator* and *Capillaria* are nematodes, *Trichophyton*
  and *Graphium* are fungi, *Trypanosoma* and *Giardia* are protists. A wrong winner means adding
  the genus to `genus_domain_override` and running again from `taxonomy_*.rds`.
- **"Records removed as homonyms"**: the removed group must be other organisms with the same genus
  name, not the same organisms classified differently. If the same organisms were removed, extend
  the "move" rule in the script.
- **"Not validly published ... kept as clinically relevant"** and **"Synonym genera with accepted
  species"**: plausible?
- **"Removed taxa that were clinically relevant"**: every removal of a taxon with prevalence < 2 needs
  an explanation (renamed, merged, truly invalid). Unexplained removals are errors.
- **"Genera that moved to another domain"**: each move changes MO codes. Accept only moves that are
  taxonomically correct.
- **"Relevant genera with a new or empty family"**: check every row. A new family is fine if it
  reflects a real reclassification; an empty family for a relevant genus is an error.
- **Duplicate names or codes, codes with repeated elements, records without a valid code**: must
  all be 0 rows. If not, find the cause.

### 3b. Sentinel organisms

Write a small R check in `<scratch>` and verify that each of these is present with the expected
properties (current name, `status`, `domain`, `family`, prefix of `mo`). Each deviation must be
explained or fixed:

| Name | Expected |
|---|---|
| *Escherichia coli*, *Klebsiella pneumoniae* | Bacteria, Enterobacteriaceae, accepted, `B_` |
| *Staphylococcus aureus* | Bacteria, Staphylococcaceae, accepted |
| *Pseudomonas aeruginosa* | Bacteria, Pseudomonadaceae |
| *Acinetobacter baumannii* | Bacteria, Moraxellaceae |
| *Enterococcus faecium* | Bacteria, Enterococcaceae |
| *Streptococcus pneumoniae* | Bacteria, Streptococcaceae |
| *Mycobacterium tuberculosis* | Bacteria, Mycobacteriaceae |
| *Clostridioides difficile* | Bacteria, Peptostreptococcaceae |
| *Tropheryma whipplei* | present (not validly published, but clinically essential) |
| *Candida albicans* | Fungi, accepted |
| *Candida auris* | synonym, with a current name in *Candidozyma* |
| *Candidozyma auris*, *Nakaseomyces glabratus* | Fungi, accepted, prevalence < 2 |
| *Aspergillus fumigatus* | Fungi, Aspergillaceae |
| *Trichophyton rubrum* | Fungi, Arthrodermataceae, `F_` code (issue #309) |
| *Mucor* (genus) | Fungi, Mucoraceae (not empty) |
| *Blastomyces* (genus), *Blastomyces dermatitidis* | Fungi, both accepted |
| *Cryptococcus neoformans*, *Pneumocystis jirovecii* | Fungi |
| *Necator americanus* | Animalia (hookworm), not Fungi |
| *Capillaria* (genus) | Animalia |
| *Ascaris lumbricoides* | Animalia, Ascarididae |
| *Echinococcus granulosus* | Animalia, Taeniidae |
| *Plasmodium falciparum* | Protozoa, `P_` code |
| *Toxoplasma gondii* | present, accepted |
| *Leishmania donovani*, *Trypanosoma cruzi*, *Giardia* | Protozoa |
| *Graphium* | only fungal species, no butterflies (e.g. no *Graphium sarpedon*) |
| *Microsporidium* | still more than 100 species |

Also verify the most used organisms in the package's own data sets: every `mo` in
`clinical_breakpoints`, `intrinsic_resistant`, `microorganisms.groups`, `microorganisms.codes` and
`example_isolates` must exist in the new `microorganisms` (`anyNA()` checks at the end of the script).

### 3c. Comparison with the previous release

Make a table for the PR with, per domain: number of records old vs new, per rank and status. Every
change of more than 5% in a domain needs an explanation. Also report:

- the number of MO codes of existing names (same full name and domain) that changed: this must be
  close to 0, every change must be explained;
- that no new taxon received an MO code that was ever used for another taxon;
- the distribution of `prevalence` by rank, old vs new;
- the 30 most relevant (prevalence <= 1.25) new genera and species, and the 30 most relevant removed
  ones.

### 3d. Integrity rules and MO code registry

The script stops before saving if any rule in `tests/testthat/helper-microorganisms.R` is broken, including
the MO code registry rules (see `data-raw/microorganisms_files/README.md`): a code from any release since
v2.0.0 may never be given to another taxon. Never weaken a rule to get past this; fix the cause. A registered
code may only get another name (e.g. a corrected spelling) if it is listed in
`data-raw/microorganisms_files/mo_code_renames.csv`, and you may only add such rows with an empty
`approved_by`: the build will then stop on it until the human approves it, so list every proposed rename under
"Needs a human". The build never removes a released taxon: missing ones are restored, except codes it retires in
`data-raw/microorganisms_files/mo_code_retirements.csv` (other organisms with the same genus name, or a genus that
is now only known in another domain). Every new row in that file needs its own checkbox in the report, and rows
marked `NEEDS REVIEW` also under "Needs a human". After the build, `tests/testthat/microorganisms_known_defects.csv` (the known defects of the
data set before this rebuild) must be emptied: delete the rows that are solved, and the file itself if none
are left. Any remaining row must be explained under "Needs a human".

### 3e. Unit tests

After the data are saved and the package is reloaded, run the tests at the end of the script and then
the full suite with `devtools::test()`. Analyse every failure: is the test outdated because of a
real taxonomic change (then update the test, and explain why in the PR) or did the build break
something (then fix the script)? Never weaken a test just to make it pass.

## Step 4: finish the package update

1. Update `TAXONOMY_VERSION` in `R/aa_globals.R` with the values from step 1.
2. Run `devtools::document()`.
3. `NEWS.md`: add one short bullet under `### Updates` (e.g. "Updated taxonomy of microorganisms to
   GBIF/COL <release>, LPSN and MycoBank of <month year>"), no full stop at the end. The version
   number in `DESCRIPTION` and on line 1 of `NEWS.md` is set by the git pre-commit hook; check after
   the first commit that it was bumped exactly once for this PR, following `CLAUDE.md`.

## Step 5: findings report for the human

Write `<review dir>/REPORT.md`. This is the permanent, git-tracked record of this update, so that the
human can review it before merging and next year's run can learn from it. Structure:

1. **Sources**: name, release or accessed date and citation of each source.
2. **Summary**: old vs new counts per domain, rank and status (the table from 3c), and the run time.
3. **Script changes**: every change to the script, why, and which records it affected.
4. **Checklist**: one checkbox per check, in this format:

   ```markdown
   - [ ] **Homonyms of relevant genera**: 98 genera in more than one domain; all winners plausible,
     except *Xus*: added to `genus_domain_override` as Animalia (nematode, 12 species), see
     [CSV](homonyms_of_relevant_genera_check_best_domain.csv)
   ```

   Every item states the check, the finding in numbers, your decision, and a link to the CSV file in
   `<review dir>` where relevant. There is one item for each review block of 3a, one for each sentinel
   organism of 3b that deviates (and one item for all sentinels that passed), one for each comparison
   of 3c, one for the integrity rules and registry (3d), one for the unit tests (3e) and one for each changed test. Leave all boxes unticked: ticking
   them is the human's job, it means "reviewed and agreed".
5. **Needs a human**: also as checkboxes, everything you could not verify or decide, each with your
   recommendation. Be concrete, e.g. "COL lists *Enterobius vermicularis* as a synonym of ...;
   recommend to keep it accepted because ...". Keep this list honest: it is the most important part
   of the report.

## Step 6: commits and pull request

Commit in logical steps, so that the human can review them separately:

1. fixes to the script (if any), with the reason in the commit message;
2. the source and intermediate files (`data-raw/taxonomy*.rds`);
3. `data/microorganisms.rda` and the other updated data sets (`fix_old_mos()` section);
4. `R/aa_globals.R`, documentation, `NEWS.md`, updated tests;
5. `<review dir>`: `REPORT.md` and all CSV files.

Push the branch and open a PR to `main` with the title "Update taxonomy of microorganisms (YYYY-MM)".
The PR description is the findings report: copy sections 1, 2 and 3 in short form, and copy the
**Checklist** and **Needs a human** sections in full, with unticked checkboxes, so that the human can
tick them one by one in the PR itself. Make the CSV links in the PR absolute links to the files on the
branch. End with a line that the PR must not be merged before all boxes are ticked.

If the human comments on the PR or asks for changes later, update the script, run again from the
appropriate checkpoint, and update both `REPORT.md` and the PR description, so that they stay identical.

## When to stop

Stop early, and report clearly what happened, only in these situations:

- source files are missing or clearly outdated (step 0);
- insufficient memory or disk space;
- the build keeps failing for reasons outside the script, or a source format changed so much that
  the parsing must be redesigned;
- sentinel checks fail and you cannot find the cause after a thorough investigation.

In the last two situations, still write `REPORT.md` with everything you found so far and the exact
point where you stopped, push the branch and open the PR **as a draft**, with the same checklist, so
that the human can continue from there.
In all other situations, finish completely and open a normal (non-draft) PR.

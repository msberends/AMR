# Data Set with Clinical Breakpoints for SIR Interpretation

Data set containing clinical breakpoints to interpret MIC and disk
diffusion to SIR values, according to international guidelines. This
dataset contain breakpoints for humans, 7 different animal groups, and
ECOFFs.

These breakpoints are currently implemented:

- For **clinical microbiology**: EUCAST 2011-2025 and CLSI 2011-2025;

- For **veterinary microbiology**: EUCAST 2021-2025 and CLSI 2019-2025;

- For **ECOFFs** (Epidemiological Cut-off Values): EUCAST 2020-2025 and
  CLSI 2022-2025.

Use [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) to transform
MICs or disks measurements to SIR values.

## Usage

``` r
clinical_breakpoints
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 40
217 observations and 14 variables:

- `guideline`  
  Name of the guideline

- `type`  
  Breakpoint type, either "ECOFF", "animal", or "human"

- `host`  
  Host of infectious agent. This is mostly useful for veterinary
  breakpoints and is either "ECOFF", "aquatic", "cats", "cattle",
  "dogs", "horse", "human", "poultry", or "swine"

- `method`  
  Testing method, either "DISK" or "MIC"

- `site`  
  Body site for which the breakpoint must be applied, e.g. "Oral" or
  "Respiratory"

- `mo`  
  Microbial ID, see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)

- `rank_index`  
  Taxonomic rank index of `mo` from 1 (subspecies/infraspecies) to 5
  (unknown microorganism)

- `ab`  
  Antimicrobial code as used by this package, EARS-Net and WHONET, see
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)

- `ref_tbl`  
  Info about where the guideline rule can be found

- `disk_dose`  
  Dose of the used disk diffusion method

- `breakpoint_S`  
  Lowest MIC value or highest number of millimetres that leads to "S"

- `breakpoint_R`  
  Highest MIC value or lowest number of millimetres that leads to "R",
  can be `NA`

- `uti`  
  A [logical](https://rdrr.io/r/base/logical.html) value
  (`TRUE`/`FALSE`) to indicate whether the rule applies to a urinary
  tract infection (UTI)

- `is_SDD`  
  A [logical](https://rdrr.io/r/base/logical.html) value
  (`TRUE`/`FALSE`) to indicate whether the intermediate range between
  "S" and "R" should be interpreted as "SDD", instead of "I". This
  currently applies to 48 breakpoints.

## Details

### Different Types of Breakpoints

Supported types of breakpoints are ECOFF, animal, and human. ECOFF
(Epidemiological cut-off) values are used in antimicrobial
susceptibility testing to differentiate between wild-type and
non-wild-type strains of bacteria or fungi.

The default is `"human"`, which can also be set with the package option
[`AMR_breakpoint_type`](https://amr-for-r.org/reference/AMR-options.md).
Use
[`as.sir(..., breakpoint_type = ...)`](https://amr-for-r.org/reference/as.sir.md)
to interpret raw data using a specific breakpoint type, e.g.
`as.sir(..., breakpoint_type = "ECOFF")` to use ECOFFs.

### Imported From WHONET

Clinical breakpoints in this package were validated through and imported
from [WHONET](https://whonet.org), a free desktop Windows application
developed and supported by the WHO Collaborating Centre for Surveillance
of Antimicrobial Resistance. More can be read on [their
website](https://whonet.org). The developers of WHONET and this `AMR`
package have been in contact about sharing their work. We highly
appreciate their great development on the WHONET software.

Our import and reproduction script can be found here:
<https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_clinical_breakpoints.R>.

### Response From CLSI and EUCAST

The CEO of CLSI and the chairman of EUCAST have endorsed the work and
public use of this `AMR` package (and consequently the use of their
breakpoints) in June 2023, when future development of distributing
clinical breakpoints was discussed in a meeting between CLSI, EUCAST,
WHO, developers of WHONET software, and developers of this `AMR`
package.

### Download Note

This `AMR` package (and the WHONET software as well) contains rather
complex internal methods to apply the guidelines. For example, some
breakpoints must be applied on certain species groups (which are in case
of this package available through the
[microorganisms.groups](https://amr-for-r.org/reference/microorganisms.groups.md)
data set). It is important that this is considered when implementing the
breakpoints for own use.

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## See also

[intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md)

## Examples

``` r
clinical_breakpoints
#> # A tibble: 40,217 × 14
#>    guideline   type  host  method site    mo            rank_index ab   ref_tbl 
#>    <chr>       <chr> <chr> <chr>  <chr>   <mo>               <dbl> <ab> <chr>   
#>  1 EUCAST 2025 human human DISK   NA      B_ACHRMB_XYLS          2 MEM  A. xylo…
#>  2 EUCAST 2025 human human MIC    NA      B_ACHRMB_XYLS          2 MEM  A. xylo…
#>  3 EUCAST 2025 human human DISK   NA      B_ACHRMB_XYLS          2 SXT  A. xylo…
#>  4 EUCAST 2025 human human MIC    NA      B_ACHRMB_XYLS          2 SXT  A. xylo…
#>  5 EUCAST 2025 human human DISK   NA      B_ACHRMB_XYLS          2 TZP  A. xylo…
#>  6 EUCAST 2025 human human MIC    NA      B_ACHRMB_XYLS          2 TZP  A. xylo…
#>  7 EUCAST 2025 human human DISK   NA      B_ACNTB                3 AMK  Acineto…
#>  8 EUCAST 2025 human human DISK   Uncomp… B_ACNTB                3 AMK  Acineto…
#>  9 EUCAST 2025 human human MIC    NA      B_ACNTB                3 AMK  Acineto…
#> 10 EUCAST 2025 human human MIC    Uncomp… B_ACNTB                3 AMK  Acineto…
#> # ℹ 40,207 more rows
#> # ℹ 5 more variables: disk_dose <chr>, breakpoint_S <dbl>, breakpoint_R <dbl>,
#> #   uti <lgl>, is_SDD <lgl>
```

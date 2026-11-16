# Data Set with Clinical Breakpoints for SIR Interpretation

Data set containing clinical breakpoints to interpret MIC and disk
diffusion to SIR values, according to international guidelines.
Currently implemented guidelines are EUCAST (2011-2023) and CLSI
(2011-2023). Use
[`as.sir()`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md) to
transform MICs or disks measurements to SIR values.

## Usage

``` r
clinical_breakpoints
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 29
747 observations and 12 variables:

- `guideline`  
  Name of the guideline

- `type`  
  Breakpoint type, either "ECOFF", "animal", or "human"

- `method`  
  Testing method, either "DISK" or "MIC"

- `site`  
  Body site for which the breakpoint must be applied, e.g. "Oral" or
  "Respiratory"

- `mo`  
  Microbial ID, see
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md)

- `rank_index`  
  Taxonomic rank index of `mo` from 1 (subspecies/infraspecies) to 5
  (unknown microorganism)

- `ab`  
  Antibiotic code as used by this package, EARS-Net and WHONET, see
  [`as.ab()`](https://amr-for-r.org/release/2.1.0/reference/as.ab.md)

- `ref_tbl`  
  Info about where the guideline rule can be found

- `disk_dose`  
  Dose of the used disk diffusion method

- `breakpoint_S`  
  Lowest MIC value or highest number of millimetres that leads to "S"

- `breakpoint_R`  
  Highest MIC value or lowest number of millimetres that leads to "R"

- `uti`  
  A [logical](https://rdrr.io/r/base/logical.html) value
  (`TRUE`/`FALSE`) to indicate whether the rule applies to a urinary
  tract infection (UTI)

## Details

### Different types of breakpoints

Supported types of breakpoints are ECOFF, animal, and human. ECOFF
(Epidemiological cut-off) values are used in antimicrobial
susceptibility testing to differentiate between wild-type and
non-wild-type strains of bacteria or fungi.

The default is `"human"`, which can also be set with the [package
option](https://amr-for-r.org/release/2.1.0/reference/AMR-options.md)
[`AMR_breakpoint_type`](https://amr-for-r.org/release/2.1.0/reference/AMR-options.md).
Use
[`as.sir(..., breakpoint_type = ...)`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md)
to interpret raw data using a specific breakpoint type, e.g.
`as.sir(..., breakpoint_type = "ECOFF")` to use ECOFFs.

### Imported from WHONET

Clinical breakpoints in this package were validated through and imported
from [WHONET](https://whonet.org), a free desktop Windows application
developed and supported by the WHO Collaborating Centre for Surveillance
of Antimicrobial Resistance. More can be read on [their
website](https://whonet.org). The developers of WHONET and this `AMR`
package have been in contact about sharing their work. We highly
appreciate their development on the WHONET software.

### Response from CLSI and EUCAST

The CEO of CLSI and the chairman of EUCAST have endorsed the work and
public use of this `AMR` package (and consequently the use of their
breakpoints) in June 2023, when future development of distributing
clinical breakpoints was discussed in a meeting between CLSI, EUCAST,
the WHO, and developers of WHONET and the `AMR` package.

### Download

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw). They
allow for machine reading EUCAST and CLSI guidelines, which is almost
impossible with the MS Excel and PDF files distributed by EUCAST and
CLSI, though initiatives have started to overcome these burdens.

**NOTE:** this `AMR` package (and the WHONET software as well) contains
internal methods to apply the guidelines, which is rather complex. For
example, some breakpoints must be applied on certain species groups
(which are in case of this package available through the
[microorganisms.groups](https://amr-for-r.org/release/2.1.0/reference/microorganisms.groups.md)
data set). It is important that this is considered when using the
breakpoints for own use.

## See also

[intrinsic_resistant](https://amr-for-r.org/release/2.1.0/reference/intrinsic_resistant.md)

## Examples

``` r
clinical_breakpoints
```

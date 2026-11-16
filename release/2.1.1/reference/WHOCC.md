# WHOCC: WHO Collaborating Centre for Drug Statistics Methodology

All antimicrobial drugs and their official names, ATC codes, ATC groups
and defined daily dose (DDD) are included in this package, using the WHO
Collaborating Centre for Drug Statistics Methodology.

## WHOCC

This package contains **all ~550 antibiotic, antimycotic and antiviral
drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC
groups and Defined Daily Dose (DDD) from the World Health Organization
Collaborating Centre for Drug Statistics Methodology (WHOCC,
<https://www.whocc.no>) and the Pharmaceuticals Community Register of
the European Commission
(<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>).

These have become the gold standard for international drug utilisation
monitoring and research.

The WHOCC is located in Oslo at the Norwegian Institute of Public Health
and funded by the Norwegian government. The European Commission is the
executive of the European Union and promotes its general interest.

**NOTE: The WHOCC copyright does not allow use for commercial purposes,
unlike any other info from this package.** See
<https://www.whocc.no/copyright_disclaimer/.>

## Examples

``` r
as.ab("meropenem")
ab_name("J01DH02")

ab_tradenames("flucloxacillin")
```


The `AMR` package for R is a powerful tool for antimicrobial resistance (AMR) analysis. It provides extensive features for handling microbial and antimicrobial data. However, for those who work primarily in Python, we now have a more intuitive option available: the [`AMR` Python package](https://pypi.org/project/AMR/).

This Python package is a wrapper around the `AMR` R package. It uses the `rpy2` package internally. Despite the need to have R installed, Python users can now easily work with AMR data directly through Python code.

# Prerequisites

This package was only tested with a [virtual environment (venv)](https://docs.python.org/3/library/venv.html). You can set up such an environment by running:

```python
# linux and macOS:
python -m venv /path/to/new/virtual/environment

# Windows:
python -m venv C:\path\to\new\virtual\environment
```

Then you can [activate the environment](https://docs.python.org/3/library/venv.html#how-venvs-work), after which the venv is ready to work with.

# Install AMR

1. Since the Python package is available on the official [Python Package Index](https://pypi.org/project/AMR/), you can just run:

    ```bash
    pip install AMR
    ```

2. Make sure you have R installed. There is **no need to install the `AMR` R package**, as it will be installed automatically.

    For Linux:

    ```bash
    # Ubuntu / Debian
    sudo apt install r-base
    # Fedora:
    sudo dnf install R
    # CentOS/RHEL
    sudo yum install R
    ```
    
    For macOS (using [Homebrew](https://brew.sh)):
    
    ```bash
    brew install r
    ```
    
    For Windows, visit the [CRAN download page](https://cran.r-project.org) to download and install R.

# Examples of Usage

## Cleaning Taxonomy

Here’s an example that demonstrates how to clean microorganism and drug names using the `AMR` Python package:

```python
import pandas as pd
import AMR

# Sample data
data = {
    "MOs": ['E. coli', 'ESCCOL', 'esco', 'Esche coli'],
    "Drug": ['Cipro', 'CIP', 'J01MA02', 'Ciproxin']
}
df = pd.DataFrame(data)

# Use AMR functions to clean microorganism and drug names
df['MO_clean'] = AMR.mo_name(df['MOs'])
df['Drug_clean'] = AMR.ab_name(df['Drug'])

# Display the results
print(df)
```

| MOs         | Drug      | MO_clean           | Drug_clean    |
|-------------|-----------|--------------------|---------------|
| E. coli     | Cipro     | Escherichia coli   | Ciprofloxacin |
| ESCCOL      | CIP       | Escherichia coli   | Ciprofloxacin |
| esco        | J01MA02   | Escherichia coli   | Ciprofloxacin |
| Esche coli  | Ciproxin  | Escherichia coli   | Ciprofloxacin |

### Explanation

* **mo_name:** This function standardises microorganism names. Here, different variations of *Escherichia coli* (such as "E. coli", "ESCCOL", "esco", and "Esche coli") are all converted into the correct, standardised form, "Escherichia coli".

* **ab_name**: Similarly, this function standardises antimicrobial names. The different representations of ciprofloxacin (e.g., "Cipro", "CIP", "J01MA02", and "Ciproxin") are all converted to the standard name, "Ciprofloxacin".

## Calculating AMR

```python
import AMR
import pandas as pd

df = AMR.example_isolates
result = AMR.resistance(df["AMX"])
print(result)
```

```
[0.59555556]
```

## Generating Antibiograms

One of the core functions of the `AMR` package is generating an antibiogram, a table that summarises the antimicrobial susceptibility of bacterial isolates. Here’s how you can generate an antibiogram from Python:

```python
result2a = AMR.antibiogram(df[["mo", "AMX", "CIP", "TZP"]])
print(result2a)
```

| Pathogen        | Amoxicillin     | Ciprofloxacin   | Piperacillin/tazobactam  |
|-----------------|-----------------|-----------------|--------------------------|
| CoNS            | 7% (10/142)     | 73% (183/252)   | 30% (10/33)              |
| E. coli         | 50% (196/392)   | 88% (399/456)   | 94% (393/416)            |
| K. pneumoniae   | 0% (0/58)       | 96% (53/55)     | 89% (47/53)              |
| P. aeruginosa   | 0% (0/30)       | 100% (30/30)    | None                     |
| P. mirabilis    | None            | 94% (34/36)     | None                     |
| S. aureus       | 6% (8/131)      | 90% (171/191)   | None                     |
| S. epidermidis  | 1% (1/91)       | 64% (87/136)    | None                     |
| S. hominis      | None            | 80% (56/70)     | None                     |
| S. pneumoniae   | 100% (112/112)  | None            | 100% (112/112)           |


```python
result2b = AMR.antibiogram(df[["mo", "AMX", "CIP", "TZP"]], mo_transform = "gramstain")
print(result2b)
```

| Pathogen       | Amoxicillin     | Ciprofloxacin    | Piperacillin/tazobactam  |
|----------------|-----------------|------------------|--------------------------|
| Gram-negative  | 36% (226/631)   | 91% (621/684)    | 88% (565/641)            |
| Gram-positive  | 43% (305/703)   | 77% (560/724)    | 86% (296/345)            |


In this example, we generate an antibiogram by selecting various antibiotics.

## Taxonomic Data Sets Now in Python!

As a Python user, you might like that the most important data sets of the `AMR` R package, `microorganisms`, `antimicrobials`, `clinical_breakpoints`, and `example_isolates`, are now available as regular Python data frames:

```python
AMR.microorganisms
```

| mo           | fullname                           | status   | kingdom  | gbif      | gbif_parent | gbif_renamed_to | prevalence |
|--------------|------------------------------------|----------|----------|-----------|-------------|-----------------|------------|
| B_GRAMN      | (unknown Gram-negatives)           | unknown  | Bacteria | None      | None        | None            | 2.0        |
| B_GRAMP      | (unknown Gram-positives)           | unknown  | Bacteria | None      | None        | None            | 2.0        |
| B_ANAER-NEG  | (unknown anaerobic Gram-negatives) | unknown  | Bacteria | None      | None        | None            | 2.0        |
| B_ANAER-POS  | (unknown anaerobic Gram-positives) | unknown  | Bacteria | None      | None        | None            | 2.0        |
| B_ANAER      | (unknown anaerobic bacteria)       | unknown  | Bacteria | None      | None        | None            | 2.0        |
| ...          | ...                                | ...      | ...      | ...       | ...         | ...             | ...        |
| B_ZYMMN_POMC | Zymomonas pomaceae                 | accepted | Bacteria | 10744418  | 3221412     | None            | 2.0        |
| B_ZYMPH      | Zymophilus                         | synonym  | Bacteria | None      | 9475166     | None            | 2.0        |
| B_ZYMPH_PCVR | Zymophilus paucivorans             | synonym  | Bacteria | None      | None        | None            | 2.0        |
| B_ZYMPH_RFFN | Zymophilus raffinosivorans         | synonym  | Bacteria | None      | None        | None            | 2.0        |
| F_ZYZYG      | Zyzygomyces                        | unknown  | Fungi    | None      | 7581        | None            | 2.0        |

```python
AMR.antimicrobials
```

| ab  | cid         | name                 | group                      | oral_ddd | oral_units | iv_ddd | iv_units |
|-----|-------------|----------------------|----------------------------|----------|------------|--------|----------|
| AMA | 4649.0      | 4-aminosalicylic acid| Antimycobacterials         | 12.00    | g          | NaN    | None     |
| ACM | 6450012.0   | Acetylmidecamycin    | Macrolides/lincosamides    | NaN      | None       | NaN    | None     |
| ASP | 49787020.0  | Acetylspiramycin     | Macrolides/lincosamides    | NaN      | None       | NaN    | None     |
| ALS | 8954.0      | Aldesulfone sodium   | Other antibacterials       | 0.33     | g          | NaN    | None     |
| AMK | 37768.0     | Amikacin             | Aminoglycosides            | NaN      | None       | 1.0    | g        |
| ... | ...         | ...                  | ...                        | ...      | ...        | ...    | ...      |
| VIR | 11979535.0  | Virginiamycine       | Other antibacterials       | NaN      | None       | NaN    | None     |
| VOR | 71616.0     | Voriconazole         | Antifungals/antimycotics   | 0.40     | g          | 0.4    | g        |
| XBR | 72144.0     | Xibornol             | Other antibacterials       | NaN      | None       | NaN    | None     |
| ZID | 77846445.0  | Zidebactam           | Other antibacterials       | NaN      | None       | NaN    | None     |
| ZFD | NaN         | Zoliflodacin         | None                       | NaN      | None       | NaN    | None     |


# Conclusion

With the `AMR` Python package, Python users can now effortlessly call R functions from the `AMR` R package. This eliminates the need for complex `rpy2` configurations and provides a clean, easy-to-use interface for antimicrobial resistance analysis. The examples provided above demonstrate how this can be applied to typical workflows, such as standardising microorganism and antimicrobial names or calculating resistance.

By just running `import AMR`, users can seamlessly integrate the robust features of the R `AMR` package into Python workflows.

Whether you're cleaning data or analysing resistance patterns, the `AMR` Python package makes it easy to work with AMR data in Python.

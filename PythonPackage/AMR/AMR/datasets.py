import os
import sys
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr, isinstalled
import pandas as pd
import importlib.metadata as metadata

# Get the path to the virtual environment
venv_path = sys.prefix

# Define R library path within the venv
r_lib_path = os.path.join(venv_path, "R_libs")
# Ensure the R library path exists
os.makedirs(r_lib_path, exist_ok=True)

# Import base and utils
base = importr('base')
utils = importr('utils')

base.options(warn = -1)

# Override R library paths globally for the session
robjects.r(f'.Library.site <- "{r_lib_path}"')  # Replace site-specific library
base._libPaths(r_lib_path)  # Override .libPaths() as well

# Get the effective library path
r_amr_lib_path = base._libPaths()[0]

# Check if the AMR package is installed in R
if not isinstalled('AMR', lib_loc=r_amr_lib_path):
    print(f"AMR: Installing latest AMR R package to {r_amr_lib_path}...", flush=True)
    utils.install_packages('AMR', repos='https://msberends.r-universe.dev', quiet=True)

# Python package version of AMR
try:
    python_amr_version = metadata.version('AMR')
except metadata.PackageNotFoundError:
    python_amr_version = ''

# R package version of AMR
r_amr_version = robjects.r(f'as.character(packageVersion("AMR", lib.loc = "{r_lib_path}"))')[0]

# Compare R and Python package versions
if r_amr_version != python_amr_version:
    try:
        print(f"AMR: Updating AMR package in {r_amr_lib_path}...", flush=True)
        utils.install_packages('AMR', repos='https://msberends.r-universe.dev', quiet=True)
    except Exception as e:
        print(f"AMR: Could not update: {e}", flush=True)

print(f"AMR: Setting up R environment and AMR datasets...", flush=True)

# Activate the automatic conversion between R and pandas DataFrames
pandas2ri.activate()

# example_isolates
example_isolates = pandas2ri.rpy2py(robjects.r('''
df <- AMR::example_isolates
df[] <- lapply(df, function(x) {
    if (inherits(x, c("Date", "POSIXt", "factor"))) {
        as.character(x)
    } else {
        x
    }
})
df <- df[, !sapply(df, is.list)]
df
'''))
example_isolates['date'] = pd.to_datetime(example_isolates['date'])

# microorganisms
microorganisms = pandas2ri.rpy2py(robjects.r('AMR::microorganisms[, !sapply(AMR::microorganisms, is.list)]'))
antimicrobials = pandas2ri.rpy2py(robjects.r('AMR::antimicrobials[, !sapply(AMR::antimicrobials, is.list)]'))
clinical_breakpoints = pandas2ri.rpy2py(robjects.r('AMR::clinical_breakpoints[, !sapply(AMR::clinical_breakpoints, is.list)]'))

base.options(warn = 0)

print(f"AMR: Done.", flush=True)

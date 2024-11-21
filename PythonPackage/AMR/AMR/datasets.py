BLUE = '\033[94m'
GREEN = '\033[32m'
RESET = '\033[0m'

print(f"{BLUE}AMR:{RESET} Setting up R environment and AMR datasets...", flush=True)


import os
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr, isinstalled
import pandas as pd
import importlib.metadata as metadata

# Get the path to the virtual environment
venv_path = os.getenv('VIRTUAL_ENV')  # Path to the active virtual environment
if not venv_path:
    raise EnvironmentError("No virtual environment detected. Please activate your Python virtual environment.")

# Define R library path within the venv
r_lib_path = os.path.join(venv_path, "R_libs")
# Ensure the R library path exists
os.makedirs(r_lib_path, exist_ok=True)
# Set the R library path in .libPaths
base = importr('base')
base._libPaths(r_lib_path)

# Check if the AMR package is installed in R
if not isinstalled('AMR'):
    utils = importr('utils')
    utils.install_packages('AMR', repos='https://msberends.r-universe.dev', quiet=True)

# Python package version of AMR
try:
    python_amr_version = metadata.version('AMR')
except metadata.PackageNotFoundError:
    python_amr_version = None

# R package version of AMR
r_amr_version = robjects.r(f'as.character(packageVersion("AMR", lib.loc = "{r_lib_path}"))')[0]

# Compare R and Python package versions
if python_amr_version and r_amr_version != python_amr_version:
    try:
        utils = importr('utils')
        utils.install_packages('AMR', repos='https://msberends.r-universe.dev', quiet=True)
    except Exception as e:
        print(f"{BLUE}AMR:{RESET} Could not update: {e}{RESET}", flush=True)

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
antibiotics = pandas2ri.rpy2py(robjects.r('AMR::antibiotics[, !sapply(AMR::antibiotics, is.list)]'))
clinical_breakpoints = pandas2ri.rpy2py(robjects.r('AMR::clinical_breakpoints[, !sapply(AMR::clinical_breakpoints, is.list)]'))

print(f"{BLUE}AMR:{RESET} {GREEN}Done.{RESET}", flush=True)

import os
import sys
import pandas as pd
import importlib.metadata as metadata

# Get the path to the virtual environment
venv_path = sys.prefix
r_lib_path = os.path.join(venv_path, "R_libs")
os.makedirs(r_lib_path, exist_ok=True)

# Set environment variable before importing rpy2
os.environ['R_LIBS_SITE'] = r_lib_path

from rpy2 import robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.packages import importr, isinstalled

# Import base and utils
base = importr('base')
utils = importr('utils')

base.options(warn=-1)

# Ensure library paths explicitly
base._libPaths(r_lib_path)

# Check if the AMR package is installed in R
if not isinstalled('AMR', lib_loc=r_lib_path):
    print(f"AMR: Installing latest AMR R package to {r_lib_path}...", flush=True)
    utils.install_packages('AMR', repos='beta.amr-for-r.org', quiet=True)

# Retrieve Python AMR version
try:
    python_amr_version = str(metadata.version('AMR'))
except metadata.PackageNotFoundError:
    python_amr_version = str('')

# Retrieve R AMR version
r_amr_version = robjects.r(f'as.character(packageVersion("AMR", lib.loc = "{r_lib_path}"))')
r_amr_version = str(r_amr_version[0])

# Compare R and Python package versions
if r_amr_version != python_amr_version:
    try:
        print(f"AMR: Updating AMR package in {r_lib_path}...", flush=True)
        utils.install_packages('AMR', repos='beta.amr-for-r.org', quiet=True)
    except Exception as e:
        print(f"AMR: Could not update: {e}", flush=True)

print(f"AMR: Setting up R environment and AMR datasets...", flush=True)

# Activate the automatic conversion between R and pandas DataFrames
with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
    # example_isolates
    example_isolates = robjects.r('''
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
    ''')
    example_isolates['date'] = pd.to_datetime(example_isolates['date'])

    # microorganisms
    microorganisms = robjects.r('AMR::microorganisms[, !sapply(AMR::microorganisms, is.list)]')
    antimicrobials = robjects.r('AMR::antimicrobials[, !sapply(AMR::antimicrobials, is.list)]')
    clinical_breakpoints = robjects.r('AMR::clinical_breakpoints[, !sapply(AMR::clinical_breakpoints, is.list)]')

base.options(warn = 0)

print(f"AMR: Done.", flush=True)

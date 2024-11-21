#!/bin/bash

# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# Clean up
rm -rf ../PythonPackage/AMR/*
mkdir -p ../PythonPackage/AMR/AMR

# Output Python file
setup_file="../PythonPackage/AMR/setup.py"
functions_file="../PythonPackage/AMR/AMR/functions.py"
datasets_file="../PythonPackage/AMR/AMR/datasets.py"
init_file="../PythonPackage/AMR/AMR/__init__.py"
description_file="../DESCRIPTION"

# Write header to the datasets Python file, including the convert_to_python function
cat <<EOL > "$datasets_file"
BLUE = '\033[94m'
GREEN = '\033[32m'
RESET = '\033[0m'

import os
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr, isinstalled
import pandas as pd
import importlib.metadata as metadata

# Get the path to the virtual environment
venv_path = os.getenv('VIRTUAL_ENV')  # Path to the active virtual environment

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
    python_amr_version = ''

# R package version of AMR
r_amr_version = robjects.r(f'as.character(packageVersion("AMR", lib.loc = "{r_lib_path}"))')[0]

# Compare R and Python package versions
if r_amr_version != python_amr_version:
    try:
        print(f"{BLUE}AMR:{RESET} Updating package version{RESET}", flush=True)
        utils = importr('utils')
        utils.install_packages('AMR', repos='https://msberends.r-universe.dev', quiet=True)
    except Exception as e:
        print(f"{BLUE}AMR:{RESET} Could not update: {e}{RESET}", flush=True)

print(f"{BLUE}AMR:{RESET} Setting up R environment and AMR datasets...", flush=True)

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
EOL

echo "from .datasets import example_isolates" >> $init_file
echo "from .datasets import microorganisms" >> $init_file
echo "from .datasets import antibiotics" >> $init_file
echo "from .datasets import clinical_breakpoints" >> $init_file


# Write header to the functions Python file, including the convert_to_python function
cat <<EOL > "$functions_file"
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, FactorVector, IntVector, FloatVector, DataFrame
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np

# Activate automatic conversion between R data frames and pandas data frames
pandas2ri.activate()

# Import the AMR R package
amr_r = importr('AMR')

def convert_to_python(r_output):
    # Check if it's a StrVector (R character vector)
    if isinstance(r_output, StrVector):
        return list(r_output)  # Convert to a Python list of strings
    
    # Check if it's a FactorVector (R factor)
    elif isinstance(r_output, FactorVector):
        return list(r_output)  # Convert to a list of integers (factor levels)
    
    # Check if it's an IntVector or FloatVector (numeric R vectors)
    elif isinstance(r_output, (IntVector, FloatVector)):
        return list(r_output)  # Convert to a Python list of integers or floats
    
    # Check if it's a pandas-compatible R data frame
    elif isinstance(r_output, pd.DataFrame):
        return r_output  # Return as pandas DataFrame (already converted by pandas2ri)
    elif isinstance(r_output, DataFrame):
        return pandas2ri.rpy2py(r_output) # Return as pandas DataFrame

    # Check if the input is a NumPy array and has a string data type
    if isinstance(r_output, np.ndarray) and np.issubdtype(r_output.dtype, np.str_):
        return r_output.tolist()  # Convert to a regular Python list
    
    # Fall-back
    return r_output
EOL

# Directory where the .Rd files are stored (update path as needed)
rd_dir="../man"

# Iterate through each .Rd file in the man directory
for rd_file in "$rd_dir"/*.Rd; do
    # Extract function names and their arguments from the .Rd files
    awk '
    BEGIN {
        usage_started = 0
    }
    
    # Detect the start of the \usage block
    /^\\usage\{/ {
        usage_started = 1
    }

    # Detect the end of the \usage block
    usage_started && /^\}/ {
        usage_started = 0
    }

    # Process lines within the \usage block that look like function calls
    usage_started && /^[a-zA-Z_]+/ {
        func_line = $0
        func_line_py = $0

        # Extract the function name (up to the first parenthesis)
        sub(/\(.*/, "", func_line)
        func_name = func_line
        func_name_py = func_name

        # Replace dots with underscores in Python function names
        gsub(/\./, "_", func_name_py)

        # Extract the arguments (inside the parentheses)
        sub(/^[^(]+\(/, "", $0)
        sub(/\).*/, "", $0)
        func_args = $0

        # Count the number of arguments
        arg_count = split(func_args, arg_array, ",")

        # Handle "..." arguments (convert them to *args, **kwargs in Python)
        gsub("\\.\\.\\.", "*args, **kwargs", func_args)

        # Remove default values from arguments
        gsub(/ = [^,]+/, "", func_args)

        # If no arguments, skip the function (dont print it)
        if (arg_count == 0) {
            func_args = "*args, **kwargs"
        }

        # If more than 1 argument, replace the 2nd to nth arguments with *args, **kwargs
        if (arg_count > 1) {
            first_arg = arg_array[1]
            func_args = first_arg ", *args, **kwargs"
        }
        if (arg_array[1] == "...") {
            func_args = "*args, **kwargs"
        }

        # Skip functions where func_name_py is identical to func_args
        if (func_name_py == func_args) {
            next
        }

        # Skip functions matching the regex pattern
        if (func_name_py ~ /^(x |facet|scale|set|get|NA_|microorganisms|antibiotics|clinical_breakpoints|example_isolates)/) {
            next
        }

        # Replace TRUE/FALSE/NULL
        gsub("TRUE", "True", func_args)
        gsub("FALSE", "False", func_args)
        gsub("NULL", "None", func_args)

        # Write the Python function definition to the output file
        print "def " func_name_py "(" func_args "):" >> "'"$functions_file"'"
        print "    \"\"\"See our website of the R package for the manual: https://msberends.github.io/AMR/index.html\"\"\"" >> "'"$functions_file"'"
        print "    return convert_to_python(amr_r." func_name_py "(" func_args "))" >> "'"$functions_file"'"
        
        print "from .functions import " func_name_py >> "'"$init_file"'"
    }
    ' "$rd_file"
done

# Output completion message
echo "Python wrapper functions generated in $functions_file."
echo "Python wrapper functions listed in $init_file."

cp ../vignettes/AMR_for_Python.Rmd ../PythonPackage/AMR/README.md
sed -i '1,/^# Introduction$/d' ../PythonPackage/AMR/README.md
echo "README copied"

# Extract the relevant fields from DESCRIPTION
version=$(grep "^Version:" "$description_file" | awk '{print $2}')

# Write the setup.py file
cat <<EOL > "$setup_file"
from setuptools import setup, find_packages

setup(
    name='AMR',
    version='$version',
    packages=find_packages(),
    install_requires=[
        'rpy2',
        'numpy',
        'pandas',
    ],
    author='Dr. Matthijs Berends',
    author_email='m.s.berends@umcg.nl',
    description='A Python wrapper for the AMR R package',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/msberends/AMR',
    project_urls={
        'Bug Tracker': 'https://github.com/msberends/AMR/issues',
    },
    license='GPL 2',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
EOL

# Output completion message
echo "setup.py has been generated in $setup_file."

cd ../PythonPackage/AMR
pip3 install build
python3 -m build
# python3 setup.py sdist bdist_wheel


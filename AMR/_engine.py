import os
import sys
import importlib.metadata as metadata

# Get the path to the virtual environment
venv_path = sys.prefix
r_lib_path = os.path.join(venv_path, "R_libs")
os.makedirs(r_lib_path, exist_ok=True)

# Set environment variable before importing rpy2
os.environ['R_LIBS_SITE'] = r_lib_path

from rpy2 import robjects
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr, isinstalled

# Import base and utils once
base = importr('base')
utils = importr('utils')

# Silence R console output entirely
robjects.r('suppressMessages(suppressWarnings(sink(tempfile())))')
base._libPaths(r_lib_path)

_installed_source = None

def _r_version():
    """Return the currently installed AMR R package version, or None."""
    try:
        return str(robjects.r(
            f'as.character(packageVersion("AMR", lib.loc = "{r_lib_path}"))')[0])
    except Exception:
        return None

def _py_version():
    """Return the Python AMR package version from metadata, or empty string."""
    try:
        return str(metadata.version('AMR'))
    except metadata.PackageNotFoundError:
        return ''

def _install_cran():
    """Install AMR from CRAN into the isolated library."""
    print("AMR: Installing from CRAN...", flush=True)
    utils.install_packages(
        'AMR',
        repos='https://cloud.r-project.org',
        lib=r_lib_path,
        quiet=True
    )

def _install_github():
    """Install AMR development version from GitHub into the isolated library."""
    print("AMR: Installing development version from GitHub...", flush=True)
    utils.install_packages(
        StrVector(['remotes', 'desc']),
        repos='https://cloud.r-project.org',
        lib=r_lib_path,
        quiet=True
    )
    remotes = importr('remotes', lib_loc=r_lib_path)
    remotes.install_github('msberends/AMR', lib=r_lib_path, quiet=True)

def ensure_amr(source="cran"):
    """Ensure AMR is installed from the requested source. Idempotent per source."""
    global _installed_source
    
    if _installed_source == source:
        return
    
    install_fn = _install_github if source == "github" else _install_cran
    
    if not isinstalled('AMR', lib_loc=r_lib_path):
        install_fn()
    else:
        # Check for version mismatch and update if needed
        r_ver = _r_version()
        py_ver = _py_version()
        if r_ver != py_ver:
            try:
                install_fn()
            except Exception as e:
                print(f"AMR: Could not update ({e})", flush=True)
    
    print(f"AMR: R package version {_r_version()} ready.", flush=True)
    _installed_source = source

def restore_sink():
    """Restore R console output after setup is complete."""
    try:
        robjects.r('sink()')
    except Exception:
        pass

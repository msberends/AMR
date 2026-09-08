import functools
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, FactorVector, IntVector, FloatVector, DataFrame
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
import pandas as pd
import numpy as np

from ._engine import ensure_amr

# Ensure AMR is available before importing it in R
ensure_amr("cran")
amr_r = importr('AMR')

def convert_to_r(value):
    """Convert Python lists/tuples to typed R vectors.

    rpy2's default_converter passes Python lists to R as R lists, not as
    character/numeric vectors. This causes element-wise type-check functions
    such as is.mic(), is.sir(), and is.disk() to return a logical vector
    rather than a single logical, breaking R's scalar && operator.

    This helper converts Python lists and tuples to the appropriate R vector
    type based on the element types, so R always receives a proper vector."""
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return StrVector([])
        # bool must be checked before int because bool is a subclass of int
        if all(isinstance(v, bool) for v in value):
            return robjects.vectors.BoolVector(value)
        if all(isinstance(v, int) for v in value):
            return IntVector(value)
        if all(isinstance(v, float) for v in value):
            return FloatVector(value)
        if all(isinstance(v, str) for v in value):
            return StrVector(value)
        # Mixed types: coerce all to string
        return StrVector([str(v) for v in value])
    return value

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
    elif isinstance(r_output, (pd.DataFrame, DataFrame)):
        return r_output  # Return as pandas DataFrame (already converted by pandas2ri)

    # Check if the input is a NumPy array and has a string data type
    if isinstance(r_output, np.ndarray) and np.issubdtype(r_output.dtype, np.str_):
        return r_output.tolist()  # Convert to a regular Python list
    
    # Fall-back
    return r_output

def r_to_python(r_func):
    """Decorator that converts Python list/tuple inputs to typed R vectors,
    runs the rpy2 function under a localconverter, and converts the output
    to a Python type."""
    @functools.wraps(r_func)
    def wrapper(*args, **kwargs):
        args = tuple(convert_to_r(a) for a in args)
        kwargs = {k: convert_to_r(v) for k, v in kwargs.items()}
        with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
            return convert_to_python(r_func(*args, **kwargs))
    return wrapper
@r_to_python
def custom_eucast_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.custom_eucast_rules(*args, **kwargs)
@r_to_python
def ab_class(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_class(*args, **kwargs)
@r_to_python
def ab_selector(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_selector(*args, **kwargs)
@r_to_python
def ab_from_text(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_from_text(*args, **kwargs)
@r_to_python
def ab_name(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_name(x, *args, **kwargs)
@r_to_python
def ab_cid(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_cid(x, *args, **kwargs)
@r_to_python
def ab_synonyms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_synonyms(x, *args, **kwargs)
@r_to_python
def ab_tradenames(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_tradenames(x, *args, **kwargs)
@r_to_python
def ab_group(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_group(x, *args, **kwargs)
@r_to_python
def ab_atc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_atc(x, *args, **kwargs)
@r_to_python
def ab_atc_group1(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_atc_group1(x, *args, **kwargs)
@r_to_python
def ab_atc_group2(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_atc_group2(x, *args, **kwargs)
@r_to_python
def ab_loinc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_loinc(x, *args, **kwargs)
@r_to_python
def ab_ddd(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_ddd(x, *args, **kwargs)
@r_to_python
def ab_ddd_units(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_ddd_units(x, *args, **kwargs)
@r_to_python
def ab_info(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_info(x, *args, **kwargs)
@r_to_python
def ab_url(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_url(x, *args, **kwargs)
@r_to_python
def ab_property(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_property(x, *args, **kwargs)
@r_to_python
def add_custom_antimicrobials(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.add_custom_antimicrobials(x)
@r_to_python
def clear_custom_antimicrobials(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.clear_custom_antimicrobials(*args, **kwargs)
@r_to_python
def add_custom_microorganisms(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.add_custom_microorganisms(x)
@r_to_python
def clear_custom_microorganisms(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.clear_custom_microorganisms(*args, **kwargs)
@r_to_python
def age(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.age(x, *args, **kwargs)
@r_to_python
def age_groups(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.age_groups(x, *args, **kwargs)
@r_to_python
def all_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_sir(*args, **kwargs)
@r_to_python
def all_sir_predictors(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_sir_predictors(*args, **kwargs)
@r_to_python
def all_mic(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_mic(*args, **kwargs)
@r_to_python
def all_mic_predictors(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_mic_predictors(*args, **kwargs)
@r_to_python
def all_disk(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_disk(*args, **kwargs)
@r_to_python
def all_disk_predictors(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_disk_predictors(*args, **kwargs)
@r_to_python
def step_mic_log2(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.step_mic_log2(*args, **kwargs)
@r_to_python
def step_sir_numeric(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.step_sir_numeric(*args, **kwargs)
@r_to_python
def amr_course(github_repo, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.amr_course(github_repo, *args, **kwargs)
@r_to_python
def wisca(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.wisca(*args, **kwargs)
@r_to_python
def antibiogram(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.antibiogram(*args, **kwargs)
@r_to_python
def retrieve_wisca_parameters(wisca_model, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.retrieve_wisca_parameters(wisca_model, *args, **kwargs)
@r_to_python
def wisca_plot(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.wisca_plot(*args, **kwargs)
@r_to_python
def aminoglycosides(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.aminoglycosides(*args, **kwargs)
@r_to_python
def aminopenicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.aminopenicillins(only_sir_columns = False, *args, **kwargs)
@r_to_python
def antifungals(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.antifungals(only_sir_columns = False, *args, **kwargs)
@r_to_python
def antimycobacterials(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.antimycobacterials(only_sir_columns = False, *args, **kwargs)
@r_to_python
def betalactams(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.betalactams(*args, **kwargs)
@r_to_python
def betalactams_with_inhibitor(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.betalactams_with_inhibitor(only_sir_columns = False, *args, **kwargs)
@r_to_python
def carbapenems(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.carbapenems(*args, **kwargs)
@r_to_python
def cephalosporins(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins(*args, **kwargs)
@r_to_python
def cephalosporins_1st(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins_1st(only_sir_columns = False, *args, **kwargs)
@r_to_python
def cephalosporins_2nd(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins_2nd(only_sir_columns = False, *args, **kwargs)
@r_to_python
def cephalosporins_3rd(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins_3rd(*args, **kwargs)
@r_to_python
def cephalosporins_4th(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins_4th(only_sir_columns = False, *args, **kwargs)
@r_to_python
def cephalosporins_5th(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.cephalosporins_5th(only_sir_columns = False, *args, **kwargs)
@r_to_python
def fluoroquinolones(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.fluoroquinolones(*args, **kwargs)
@r_to_python
def glycopeptides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.glycopeptides(only_sir_columns = False, *args, **kwargs)
@r_to_python
def ionophores(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ionophores(only_sir_columns = False, *args, **kwargs)
@r_to_python
def isoxazolylpenicillins(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.isoxazolylpenicillins(*args, **kwargs)
@r_to_python
def lincosamides(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.lincosamides(*args, **kwargs)
@r_to_python
def lipoglycopeptides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.lipoglycopeptides(only_sir_columns = False, *args, **kwargs)
@r_to_python
def macrolides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.macrolides(only_sir_columns = False, *args, **kwargs)
@r_to_python
def monobactams(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.monobactams(only_sir_columns = False, *args, **kwargs)
@r_to_python
def nitrofurans(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.nitrofurans(only_sir_columns = False, *args, **kwargs)
@r_to_python
def oxazolidinones(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.oxazolidinones(only_sir_columns = False, *args, **kwargs)
@r_to_python
def penicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.penicillins(only_sir_columns = False, *args, **kwargs)
@r_to_python
def peptides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.peptides(only_sir_columns = False, *args, **kwargs)
@r_to_python
def phenicols(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.phenicols(only_sir_columns = False, *args, **kwargs)
@r_to_python
def phosphonics(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.phosphonics(only_sir_columns = False, *args, **kwargs)
@r_to_python
def polymyxins(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.polymyxins(*args, **kwargs)
@r_to_python
def quinolones(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.quinolones(*args, **kwargs)
@r_to_python
def rifamycins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.rifamycins(only_sir_columns = False, *args, **kwargs)
@r_to_python
def spiropyrimidinetriones(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.spiropyrimidinetriones(only_sir_columns = False, *args, **kwargs)
@r_to_python
def streptogramins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.streptogramins(only_sir_columns = False, *args, **kwargs)
@r_to_python
def sulfonamides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.sulfonamides(only_sir_columns = False, *args, **kwargs)
@r_to_python
def tetracyclines(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.tetracyclines(*args, **kwargs)
@r_to_python
def trimethoprims(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.trimethoprims(only_sir_columns = False, *args, **kwargs)
@r_to_python
def ureidopenicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ureidopenicillins(only_sir_columns = False, *args, **kwargs)
@r_to_python
def amr_class(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.amr_class(*args, **kwargs)
@r_to_python
def amr_selector(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.amr_selector(*args, **kwargs)
@r_to_python
def administrable_per_os(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.administrable_per_os(only_sir_columns = False, *args, **kwargs)
@r_to_python
def administrable_iv(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.administrable_iv(only_sir_columns = False, *args, **kwargs)
@r_to_python
def not_intrinsic_resistant(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.not_intrinsic_resistant(*args, **kwargs)
@r_to_python
def as_ab(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_ab(*args, **kwargs)
@r_to_python
def is_ab(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_ab(x)
@r_to_python
def ab_reset_session(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ab_reset_session(*args, **kwargs)
@r_to_python
def as_av(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_av(x, *args, **kwargs)
@r_to_python
def is_av(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_av(x)
@r_to_python
def as_disk(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_disk(x, *args, **kwargs)
@r_to_python
def is_disk(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_disk(x)
@r_to_python
def as_mic(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_mic(x, *args, **kwargs)
@r_to_python
def is_mic(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_mic(x)
@r_to_python
def rescale_mic(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.rescale_mic(*args, **kwargs)
@r_to_python
def mic_p50(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mic_p50(x, *args, **kwargs)
@r_to_python
def mic_p90(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mic_p90(x, *args, **kwargs)
@r_to_python
def as_mo(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_mo(*args, **kwargs)
@r_to_python
def is_mo(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_mo(x)
@r_to_python
def mo_uncertainties(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_uncertainties(*args, **kwargs)
@r_to_python
def mo_renamed(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_renamed(*args, **kwargs)
@r_to_python
def mo_failures(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_failures(*args, **kwargs)
@r_to_python
def mo_reset_session(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_reset_session(*args, **kwargs)
@r_to_python
def mo_cleaning_regex(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_cleaning_regex(*args, **kwargs)
@r_to_python
def as_sir(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.as_sir(x, *args, **kwargs)
@r_to_python
def is_sir(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_sir(x)
@r_to_python
def is_sir_eligible(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_sir_eligible(x, *args, **kwargs)
@r_to_python
def sir_interpretation_history(clean):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.sir_interpretation_history(clean)
@r_to_python
def atc_online_property(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.atc_online_property(*args, **kwargs)
@r_to_python
def atc_online_groups(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.atc_online_groups(atc_code, *args, **kwargs)
@r_to_python
def atc_online_ddd(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.atc_online_ddd(atc_code, *args, **kwargs)
@r_to_python
def atc_online_ddd_units(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.atc_online_ddd_units(atc_code, *args, **kwargs)
@r_to_python
def av_from_text(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_from_text(*args, **kwargs)
@r_to_python
def av_name(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_name(x, *args, **kwargs)
@r_to_python
def av_cid(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_cid(x, *args, **kwargs)
@r_to_python
def av_synonyms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_synonyms(x, *args, **kwargs)
@r_to_python
def av_tradenames(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_tradenames(x, *args, **kwargs)
@r_to_python
def av_group(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_group(x, *args, **kwargs)
@r_to_python
def av_atc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_atc(x, *args, **kwargs)
@r_to_python
def av_loinc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_loinc(x, *args, **kwargs)
@r_to_python
def av_ddd(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_ddd(x, *args, **kwargs)
@r_to_python
def av_ddd_units(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_ddd_units(x, *args, **kwargs)
@r_to_python
def av_info(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_info(x, *args, **kwargs)
@r_to_python
def av_url(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_url(x, *args, **kwargs)
@r_to_python
def av_property(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.av_property(x, *args, **kwargs)
@r_to_python
def availability(tbl, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.availability(tbl, *args, **kwargs)
@r_to_python
def bug_drug_combinations(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.bug_drug_combinations(*args, **kwargs)
@r_to_python
def count_resistant(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_resistant(*args, **kwargs)
@r_to_python
def count_susceptible(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_susceptible(*args, **kwargs)
@r_to_python
def count_S(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_S(*args, **kwargs)
@r_to_python
def count_SI(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_SI(*args, **kwargs)
@r_to_python
def count_I(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_I(*args, **kwargs)
@r_to_python
def count_IR(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_IR(*args, **kwargs)
@r_to_python
def count_R(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_R(*args, **kwargs)
@r_to_python
def count_all(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_all(*args, **kwargs)
@r_to_python
def n_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.n_sir(*args, **kwargs)
@r_to_python
def count_df(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.count_df(*args, **kwargs)
@r_to_python
def custom_interpretive_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.custom_interpretive_rules(*args, **kwargs)
@r_to_python
def custom_mdro_guideline(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.custom_mdro_guideline(*args, **kwargs)
@r_to_python
def export_ncbi_biosample(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.export_ncbi_biosample(*args, **kwargs)
@r_to_python
def first_isolate(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.first_isolate(*args, **kwargs)
@r_to_python
def filter_first_isolate(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.filter_first_isolate(*args, **kwargs)
@r_to_python
def g_test(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.g_test(x, *args, **kwargs)
@r_to_python
def is_new_episode(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.is_new_episode(x, *args, **kwargs)
@r_to_python
def ggplot_pca(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ggplot_pca(*args, **kwargs)
@r_to_python
def ggplot_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ggplot_sir(*args, **kwargs)
@r_to_python
def geom_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.geom_sir(*args, **kwargs)
@r_to_python
def guess_ab_col(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.guess_ab_col(*args, **kwargs)
@r_to_python
def interpretive_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.interpretive_rules(*args, **kwargs)
@r_to_python
def eucast_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.eucast_rules(*args, **kwargs)
@r_to_python
def clsi_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.clsi_rules(*args, **kwargs)
@r_to_python
def eucast_dosage(ab, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.eucast_dosage(ab, *args, **kwargs)
@r_to_python
def italicise_taxonomy(string, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.italicise_taxonomy(string, *args, **kwargs)
@r_to_python
def italicize_taxonomy(string, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.italicize_taxonomy(string, *args, **kwargs)
@r_to_python
def inner_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.inner_join_microorganisms(x, *args, **kwargs)
@r_to_python
def left_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.left_join_microorganisms(x, *args, **kwargs)
@r_to_python
def right_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.right_join_microorganisms(x, *args, **kwargs)
@r_to_python
def full_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.full_join_microorganisms(x, *args, **kwargs)
@r_to_python
def semi_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.semi_join_microorganisms(x, *args, **kwargs)
@r_to_python
def anti_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.anti_join_microorganisms(x, *args, **kwargs)
@r_to_python
def key_antimicrobials(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.key_antimicrobials(*args, **kwargs)
@r_to_python
def all_antimicrobials(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.all_antimicrobials(x = None, *args, **kwargs)
@r_to_python
def kurtosis(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.kurtosis(x, *args, **kwargs)
@r_to_python
def like(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.like(x, *args, **kwargs)
@r_to_python
def mdro(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mdro(*args, **kwargs)
@r_to_python
def brmo(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.brmo(x = None, *args, **kwargs)
@r_to_python
def mrgn(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mrgn(x = None, *args, **kwargs)
@r_to_python
def mdr_tb(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mdr_tb(x = None, *args, **kwargs)
@r_to_python
def mdr_cmi2012(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mdr_cmi2012(x = None, *args, **kwargs)
@r_to_python
def eucast_exceptional_phenotypes(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.eucast_exceptional_phenotypes(*args, **kwargs)
@r_to_python
def mean_amr_distance(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mean_amr_distance(x, *args, **kwargs)
@r_to_python
def amr_distance_from_row(amr_distance, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.amr_distance_from_row(amr_distance, *args, **kwargs)
@r_to_python
def mo_matching_score(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_matching_score(x, *args, **kwargs)
@r_to_python
def mo_name(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_name(*args, **kwargs)
@r_to_python
def mo_fullname(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_fullname(*args, **kwargs)
@r_to_python
def mo_shortname(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_shortname(*args, **kwargs)
@r_to_python
def mo_subspecies(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_subspecies(*args, **kwargs)
@r_to_python
def mo_species(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_species(*args, **kwargs)
@r_to_python
def mo_genus(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_genus(*args, **kwargs)
@r_to_python
def mo_family(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_family(*args, **kwargs)
@r_to_python
def mo_order(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_order(*args, **kwargs)
@r_to_python
def mo_class(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_class(*args, **kwargs)
@r_to_python
def mo_phylum(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_phylum(*args, **kwargs)
@r_to_python
def mo_kingdom(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_kingdom(*args, **kwargs)
@r_to_python
def mo_domain(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_domain(*args, **kwargs)
@r_to_python
def mo_type(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_type(*args, **kwargs)
@r_to_python
def mo_status(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_status(*args, **kwargs)
@r_to_python
def mo_pathogenicity(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_pathogenicity(*args, **kwargs)
@r_to_python
def mo_gramstain(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_gramstain(*args, **kwargs)
@r_to_python
def mo_is_gram_negative(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_is_gram_negative(*args, **kwargs)
@r_to_python
def mo_is_gram_positive(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_is_gram_positive(*args, **kwargs)
@r_to_python
def mo_is_yeast(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_is_yeast(*args, **kwargs)
@r_to_python
def mo_is_intrinsic_resistant(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_is_intrinsic_resistant(*args, **kwargs)
@r_to_python
def mo_oxygen_tolerance(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_oxygen_tolerance(*args, **kwargs)
@r_to_python
def mo_is_anaerobic(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_is_anaerobic(*args, **kwargs)
@r_to_python
def mo_morphology(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_morphology(*args, **kwargs)
@r_to_python
def mo_snomed(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_snomed(*args, **kwargs)
@r_to_python
def mo_ref(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_ref(*args, **kwargs)
@r_to_python
def mo_authors(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_authors(*args, **kwargs)
@r_to_python
def mo_year(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_year(*args, **kwargs)
@r_to_python
def mo_lpsn(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_lpsn(*args, **kwargs)
@r_to_python
def mo_mycobank(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_mycobank(*args, **kwargs)
@r_to_python
def mo_gbif(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_gbif(*args, **kwargs)
@r_to_python
def mo_rank(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_rank(*args, **kwargs)
@r_to_python
def mo_taxonomy(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_taxonomy(*args, **kwargs)
@r_to_python
def mo_synonyms(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_synonyms(*args, **kwargs)
@r_to_python
def mo_current(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_current(x, *args, **kwargs)
@r_to_python
def mo_group_members(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_group_members(*args, **kwargs)
@r_to_python
def mo_info(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_info(*args, **kwargs)
@r_to_python
def mo_url(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_url(*args, **kwargs)
@r_to_python
def mo_property(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.mo_property(*args, **kwargs)
@r_to_python
def pca(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.pca(*args, **kwargs)
@r_to_python
def theme_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.theme_sir(*args, **kwargs)
@r_to_python
def labels_sir_count(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.labels_sir_count(*args, **kwargs)
@r_to_python
def resistance(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.resistance(*args, **kwargs)
@r_to_python
def susceptibility(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.susceptibility(*args, **kwargs)
@r_to_python
def sir_confidence_interval(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.sir_confidence_interval(*args, **kwargs)
@r_to_python
def proportion_R(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_R(*args, **kwargs)
@r_to_python
def proportion_IR(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_IR(*args, **kwargs)
@r_to_python
def proportion_I(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_I(*args, **kwargs)
@r_to_python
def proportion_SI(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_SI(*args, **kwargs)
@r_to_python
def proportion_S(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_S(*args, **kwargs)
@r_to_python
def proportion_df(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.proportion_df(*args, **kwargs)
@r_to_python
def sir_df(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.sir_df(*args, **kwargs)
@r_to_python
def random_mic(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.random_mic(*args, **kwargs)
@r_to_python
def random_disk(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.random_disk(*args, **kwargs)
@r_to_python
def random_sir(size = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.random_sir(size = None, *args, **kwargs)
@r_to_python
def resistance_predict(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.resistance_predict(*args, **kwargs)
@r_to_python
def sir_predict(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.sir_predict(*args, **kwargs)
@r_to_python
def ggplot_sir_predict(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.ggplot_sir_predict(*args, **kwargs)
@r_to_python
def skewness(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.skewness(x, *args, **kwargs)
@r_to_python
def top_n_microorganisms(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.top_n_microorganisms(*args, **kwargs)
@r_to_python
def reset_AMR_locale(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.reset_AMR_locale(*args, **kwargs)
@r_to_python
def translate_AMR(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return amr_r.translate_AMR(x, *args, **kwargs)

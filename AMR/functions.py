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
def ab_class(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_class(*args, **kwargs))
def ab_selector(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_selector(*args, **kwargs))
def ab_from_text(text, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_from_text(text, *args, **kwargs))
def ab_name(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_name(x, *args, **kwargs))
def ab_cid(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_cid(x, *args, **kwargs))
def ab_synonyms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_synonyms(x, *args, **kwargs))
def ab_tradenames(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_tradenames(x, *args, **kwargs))
def ab_group(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_group(x, *args, **kwargs))
def ab_atc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_atc(x, *args, **kwargs))
def ab_atc_group1(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_atc_group1(x, *args, **kwargs))
def ab_atc_group2(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_atc_group2(x, *args, **kwargs))
def ab_loinc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_loinc(x, *args, **kwargs))
def ab_ddd(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_ddd(x, *args, **kwargs))
def ab_ddd_units(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_ddd_units(x, *args, **kwargs))
def ab_info(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_info(x, *args, **kwargs))
def ab_url(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_url(x, *args, **kwargs))
def ab_property(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_property(x, *args, **kwargs))
def add_custom_antimicrobials(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.add_custom_antimicrobials(x))
def clear_custom_antimicrobials(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.clear_custom_antimicrobials(*args, **kwargs))
def add_custom_microorganisms(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.add_custom_microorganisms(x))
def clear_custom_microorganisms(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.clear_custom_microorganisms(*args, **kwargs))
def age(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.age(x, *args, **kwargs))
def age_groups(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.age_groups(x, *args, **kwargs))
def antibiogram(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.antibiogram(x, *args, **kwargs))
def wisca(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.wisca(x, *args, **kwargs))
def retrieve_wisca_parameters(wisca_model, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.retrieve_wisca_parameters(wisca_model, *args, **kwargs))
def aminoglycosides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.aminoglycosides(only_sir_columns = False, *args, **kwargs))
def aminopenicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.aminopenicillins(only_sir_columns = False, *args, **kwargs))
def antifungals(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.antifungals(only_sir_columns = False, *args, **kwargs))
def antimycobacterials(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.antimycobacterials(only_sir_columns = False, *args, **kwargs))
def betalactams(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.betalactams(only_sir_columns = False, *args, **kwargs))
def betalactams_with_inhibitor(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.betalactams_with_inhibitor(only_sir_columns = False, *args, **kwargs))
def carbapenems(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.carbapenems(only_sir_columns = False, *args, **kwargs))
def cephalosporins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins(only_sir_columns = False, *args, **kwargs))
def cephalosporins_1st(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins_1st(only_sir_columns = False, *args, **kwargs))
def cephalosporins_2nd(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins_2nd(only_sir_columns = False, *args, **kwargs))
def cephalosporins_3rd(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins_3rd(only_sir_columns = False, *args, **kwargs))
def cephalosporins_4th(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins_4th(only_sir_columns = False, *args, **kwargs))
def cephalosporins_5th(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.cephalosporins_5th(only_sir_columns = False, *args, **kwargs))
def fluoroquinolones(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.fluoroquinolones(only_sir_columns = False, *args, **kwargs))
def glycopeptides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.glycopeptides(only_sir_columns = False, *args, **kwargs))
def isoxazolylpenicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.isoxazolylpenicillins(only_sir_columns = False, *args, **kwargs))
def lincosamides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.lincosamides(only_sir_columns = False, *args, **kwargs))
def lipoglycopeptides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.lipoglycopeptides(only_sir_columns = False, *args, **kwargs))
def macrolides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.macrolides(only_sir_columns = False, *args, **kwargs))
def monobactams(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.monobactams(only_sir_columns = False, *args, **kwargs))
def nitrofurans(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.nitrofurans(only_sir_columns = False, *args, **kwargs))
def oxazolidinones(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.oxazolidinones(only_sir_columns = False, *args, **kwargs))
def penicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.penicillins(only_sir_columns = False, *args, **kwargs))
def phenicols(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.phenicols(only_sir_columns = False, *args, **kwargs))
def polymyxins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.polymyxins(only_sir_columns = False, *args, **kwargs))
def quinolones(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.quinolones(only_sir_columns = False, *args, **kwargs))
def rifamycins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.rifamycins(only_sir_columns = False, *args, **kwargs))
def streptogramins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.streptogramins(only_sir_columns = False, *args, **kwargs))
def sulfonamides(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.sulfonamides(only_sir_columns = False, *args, **kwargs))
def tetracyclines(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.tetracyclines(only_sir_columns = False, *args, **kwargs))
def trimethoprims(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.trimethoprims(only_sir_columns = False, *args, **kwargs))
def ureidopenicillins(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ureidopenicillins(only_sir_columns = False, *args, **kwargs))
def amr_class(amr_class, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.amr_class(amr_class, *args, **kwargs))
def amr_selector(filter, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.amr_selector(filter, *args, **kwargs))
def administrable_per_os(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.administrable_per_os(only_sir_columns = False, *args, **kwargs))
def administrable_iv(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.administrable_iv(only_sir_columns = False, *args, **kwargs))
def not_intrinsic_resistant(only_sir_columns = False, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.not_intrinsic_resistant(only_sir_columns = False, *args, **kwargs))
def as_ab(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_ab(x, *args, **kwargs))
def is_ab(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_ab(x))
def ab_reset_session(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ab_reset_session(*args, **kwargs))
def as_av(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_av(x, *args, **kwargs))
def is_av(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_av(x))
def as_disk(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_disk(x, *args, **kwargs))
def is_disk(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_disk(x))
def as_mic(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_mic(x, *args, **kwargs))
def is_mic(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_mic(x))
def rescale_mic(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.rescale_mic(x, *args, **kwargs))
def mic_p50(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mic_p50(x, *args, **kwargs))
def mic_p90(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mic_p90(x, *args, **kwargs))
def as_mo(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_mo(x, *args, **kwargs))
def is_mo(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_mo(x))
def mo_uncertainties(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_uncertainties(*args, **kwargs))
def mo_renamed(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_renamed(*args, **kwargs))
def mo_failures(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_failures(*args, **kwargs))
def mo_reset_session(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_reset_session(*args, **kwargs))
def mo_cleaning_regex(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_cleaning_regex(*args, **kwargs))
def as_sir(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.as_sir(x, *args, **kwargs))
def is_sir(x):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_sir(x))
def is_sir_eligible(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_sir_eligible(x, *args, **kwargs))
def sir_interpretation_history(clean):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.sir_interpretation_history(clean))
def atc_online_property(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.atc_online_property(atc_code, *args, **kwargs))
def atc_online_groups(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.atc_online_groups(atc_code, *args, **kwargs))
def atc_online_ddd(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.atc_online_ddd(atc_code, *args, **kwargs))
def atc_online_ddd_units(atc_code, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.atc_online_ddd_units(atc_code, *args, **kwargs))
def av_from_text(text, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_from_text(text, *args, **kwargs))
def av_name(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_name(x, *args, **kwargs))
def av_cid(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_cid(x, *args, **kwargs))
def av_synonyms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_synonyms(x, *args, **kwargs))
def av_tradenames(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_tradenames(x, *args, **kwargs))
def av_group(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_group(x, *args, **kwargs))
def av_atc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_atc(x, *args, **kwargs))
def av_loinc(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_loinc(x, *args, **kwargs))
def av_ddd(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_ddd(x, *args, **kwargs))
def av_ddd_units(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_ddd_units(x, *args, **kwargs))
def av_info(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_info(x, *args, **kwargs))
def av_url(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_url(x, *args, **kwargs))
def av_property(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.av_property(x, *args, **kwargs))
def availability(tbl, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.availability(tbl, *args, **kwargs))
def bug_drug_combinations(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.bug_drug_combinations(x, *args, **kwargs))
def count_resistant(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_resistant(*args, **kwargs))
def count_susceptible(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_susceptible(*args, **kwargs))
def count_S(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_S(*args, **kwargs))
def count_SI(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_SI(*args, **kwargs))
def count_I(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_I(*args, **kwargs))
def count_IR(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_IR(*args, **kwargs))
def count_R(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_R(*args, **kwargs))
def count_all(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_all(*args, **kwargs))
def n_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.n_sir(*args, **kwargs))
def count_df(data, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.count_df(data, *args, **kwargs))
def custom_eucast_rules(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.custom_eucast_rules(*args, **kwargs))
def eucast_rules(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.eucast_rules(x, *args, **kwargs))
def eucast_dosage(ab, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.eucast_dosage(ab, *args, **kwargs))
def export_ncbi_biosample(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.export_ncbi_biosample(x, *args, **kwargs))
def first_isolate(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.first_isolate(x = None, *args, **kwargs))
def filter_first_isolate(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.filter_first_isolate(x = None, *args, **kwargs))
def g_test(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.g_test(x, *args, **kwargs))
def is_new_episode(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.is_new_episode(x, *args, **kwargs))
def ggplot_pca(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ggplot_pca(x, *args, **kwargs))
def ggplot_sir(data, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ggplot_sir(data, *args, **kwargs))
def geom_sir(position = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.geom_sir(position = None, *args, **kwargs))
def guess_ab_col(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.guess_ab_col(x = None, *args, **kwargs))
def italicise_taxonomy(string, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.italicise_taxonomy(string, *args, **kwargs))
def italicize_taxonomy(string, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.italicize_taxonomy(string, *args, **kwargs))
def inner_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.inner_join_microorganisms(x, *args, **kwargs))
def left_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.left_join_microorganisms(x, *args, **kwargs))
def right_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.right_join_microorganisms(x, *args, **kwargs))
def full_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.full_join_microorganisms(x, *args, **kwargs))
def semi_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.semi_join_microorganisms(x, *args, **kwargs))
def anti_join_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.anti_join_microorganisms(x, *args, **kwargs))
def key_antimicrobials(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.key_antimicrobials(x = None, *args, **kwargs))
def all_antimicrobials(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.all_antimicrobials(x = None, *args, **kwargs))
def kurtosis(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.kurtosis(x, *args, **kwargs))
def like(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.like(x, *args, **kwargs))
def mdro(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mdro(x = None, *args, **kwargs))
def custom_mdro_guideline(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.custom_mdro_guideline(*args, **kwargs))
def brmo(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.brmo(x = None, *args, **kwargs))
def mrgn(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mrgn(x = None, *args, **kwargs))
def mdr_tb(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mdr_tb(x = None, *args, **kwargs))
def mdr_cmi2012(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mdr_cmi2012(x = None, *args, **kwargs))
def eucast_exceptional_phenotypes(x = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.eucast_exceptional_phenotypes(x = None, *args, **kwargs))
def mean_amr_distance(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mean_amr_distance(x, *args, **kwargs))
def amr_distance_from_row(amr_distance, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.amr_distance_from_row(amr_distance, *args, **kwargs))
def mo_matching_score(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_matching_score(x, *args, **kwargs))
def mo_name(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_name(x, *args, **kwargs))
def mo_fullname(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_fullname(x, *args, **kwargs))
def mo_shortname(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_shortname(x, *args, **kwargs))
def mo_subspecies(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_subspecies(x, *args, **kwargs))
def mo_species(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_species(x, *args, **kwargs))
def mo_genus(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_genus(x, *args, **kwargs))
def mo_family(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_family(x, *args, **kwargs))
def mo_order(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_order(x, *args, **kwargs))
def mo_class(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_class(x, *args, **kwargs))
def mo_phylum(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_phylum(x, *args, **kwargs))
def mo_kingdom(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_kingdom(x, *args, **kwargs))
def mo_domain(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_domain(x, *args, **kwargs))
def mo_type(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_type(x, *args, **kwargs))
def mo_status(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_status(x, *args, **kwargs))
def mo_pathogenicity(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_pathogenicity(x, *args, **kwargs))
def mo_gramstain(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_gramstain(x, *args, **kwargs))
def mo_is_gram_negative(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_is_gram_negative(x, *args, **kwargs))
def mo_is_gram_positive(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_is_gram_positive(x, *args, **kwargs))
def mo_is_yeast(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_is_yeast(x, *args, **kwargs))
def mo_is_intrinsic_resistant(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_is_intrinsic_resistant(x, *args, **kwargs))
def mo_oxygen_tolerance(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_oxygen_tolerance(x, *args, **kwargs))
def mo_is_anaerobic(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_is_anaerobic(x, *args, **kwargs))
def mo_snomed(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_snomed(x, *args, **kwargs))
def mo_ref(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_ref(x, *args, **kwargs))
def mo_authors(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_authors(x, *args, **kwargs))
def mo_year(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_year(x, *args, **kwargs))
def mo_lpsn(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_lpsn(x, *args, **kwargs))
def mo_mycobank(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_mycobank(x, *args, **kwargs))
def mo_gbif(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_gbif(x, *args, **kwargs))
def mo_rank(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_rank(x, *args, **kwargs))
def mo_taxonomy(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_taxonomy(x, *args, **kwargs))
def mo_synonyms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_synonyms(x, *args, **kwargs))
def mo_current(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_current(x, *args, **kwargs))
def mo_group_members(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_group_members(x, *args, **kwargs))
def mo_info(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_info(x, *args, **kwargs))
def mo_url(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_url(x, *args, **kwargs))
def mo_property(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.mo_property(x, *args, **kwargs))
def pca(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.pca(x, *args, **kwargs))
def theme_sir(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.theme_sir(*args, **kwargs))
def labels_sir_count(position = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.labels_sir_count(position = None, *args, **kwargs))
def resistance(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.resistance(*args, **kwargs))
def susceptibility(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.susceptibility(*args, **kwargs))
def sir_confidence_interval(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.sir_confidence_interval(*args, **kwargs))
def proportion_R(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_R(*args, **kwargs))
def proportion_IR(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_IR(*args, **kwargs))
def proportion_I(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_I(*args, **kwargs))
def proportion_SI(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_SI(*args, **kwargs))
def proportion_S(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_S(*args, **kwargs))
def proportion_df(data, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.proportion_df(data, *args, **kwargs))
def sir_df(data, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.sir_df(data, *args, **kwargs))
def random_mic(size = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.random_mic(size = None, *args, **kwargs))
def random_disk(size = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.random_disk(size = None, *args, **kwargs))
def random_sir(size = None, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.random_sir(size = None, *args, **kwargs))
def resistance_predict(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.resistance_predict(x, *args, **kwargs))
def sir_predict(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.sir_predict(x, *args, **kwargs))
def ggplot_sir_predict(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.ggplot_sir_predict(x, *args, **kwargs))
def skewness(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.skewness(x, *args, **kwargs))
def top_n_microorganisms(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.top_n_microorganisms(x, *args, **kwargs))
def reset_AMR_locale(*args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.reset_AMR_locale(*args, **kwargs))
def translate_AMR(x, *args, **kwargs):
    """Please see our website of the R package for the full manual: https://amr-for-r.org"""
    return convert_to_python(amr_r.translate_AMR(x, *args, **kwargs))

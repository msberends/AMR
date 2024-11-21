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
def ab_from_text(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_from_text(*args, **kwargs))
def ab_name(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_name(x, *args, **kwargs))
def ab_cid(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_cid(x, *args, **kwargs))
def ab_synonyms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_synonyms(x, *args, **kwargs))
def ab_tradenames(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_tradenames(x, *args, **kwargs))
def ab_group(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_group(x, *args, **kwargs))
def ab_atc(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_atc(x, *args, **kwargs))
def ab_atc_group1(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_atc_group1(x, *args, **kwargs))
def ab_atc_group2(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_atc_group2(x, *args, **kwargs))
def ab_loinc(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_loinc(x, *args, **kwargs))
def ab_ddd(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_ddd(x, *args, **kwargs))
def ab_ddd_units(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_ddd_units(x, *args, **kwargs))
def ab_info(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_info(x, *args, **kwargs))
def ab_url(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_url(x, *args, **kwargs))
def ab_property(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_property(x, *args, **kwargs))
def add_custom_antimicrobials(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.add_custom_antimicrobials(x))
def clear_custom_antimicrobials(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.clear_custom_antimicrobials(*args, **kwargs))
def add_custom_microorganisms(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.add_custom_microorganisms(x))
def clear_custom_microorganisms(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.clear_custom_microorganisms(*args, **kwargs))
def age(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.age(x, *args, **kwargs))
def age_groups(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.age_groups(x, *args, **kwargs))
def antibiogram(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.antibiogram(*args, **kwargs))
def ab_class(ab_class, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_class(ab_class, *args, **kwargs))
def ab_selector(filter, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ab_selector(filter, *args, **kwargs))
def aminoglycosides(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.aminoglycosides(only_sir_columns = False, *args, **kwargs))
def aminopenicillins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.aminopenicillins(only_sir_columns = False, *args, **kwargs))
def antifungals(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.antifungals(only_sir_columns = False, *args, **kwargs))
def antimycobacterials(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.antimycobacterials(only_sir_columns = False, *args, **kwargs))
def betalactams(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.betalactams(only_sir_columns = False, *args, **kwargs))
def carbapenems(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.carbapenems(only_sir_columns = False, *args, **kwargs))
def cephalosporins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins(only_sir_columns = False, *args, **kwargs))
def cephalosporins_1st(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins_1st(only_sir_columns = False, *args, **kwargs))
def cephalosporins_2nd(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins_2nd(only_sir_columns = False, *args, **kwargs))
def cephalosporins_3rd(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins_3rd(only_sir_columns = False, *args, **kwargs))
def cephalosporins_4th(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins_4th(only_sir_columns = False, *args, **kwargs))
def cephalosporins_5th(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.cephalosporins_5th(only_sir_columns = False, *args, **kwargs))
def fluoroquinolones(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.fluoroquinolones(only_sir_columns = False, *args, **kwargs))
def glycopeptides(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.glycopeptides(only_sir_columns = False, *args, **kwargs))
def lincosamides(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.lincosamides(only_sir_columns = False, *args, **kwargs))
def lipoglycopeptides(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.lipoglycopeptides(only_sir_columns = False, *args, **kwargs))
def macrolides(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.macrolides(only_sir_columns = False, *args, **kwargs))
def nitrofurans(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.nitrofurans(only_sir_columns = False, *args, **kwargs))
def oxazolidinones(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.oxazolidinones(only_sir_columns = False, *args, **kwargs))
def penicillins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.penicillins(only_sir_columns = False, *args, **kwargs))
def polymyxins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.polymyxins(only_sir_columns = False, *args, **kwargs))
def quinolones(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.quinolones(only_sir_columns = False, *args, **kwargs))
def rifamycins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.rifamycins(only_sir_columns = False, *args, **kwargs))
def streptogramins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.streptogramins(only_sir_columns = False, *args, **kwargs))
def tetracyclines(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.tetracyclines(only_sir_columns = False, *args, **kwargs))
def trimethoprims(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.trimethoprims(only_sir_columns = False, *args, **kwargs))
def ureidopenicillins(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ureidopenicillins(only_sir_columns = False, *args, **kwargs))
def administrable_per_os(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.administrable_per_os(only_sir_columns = False, *args, **kwargs))
def administrable_iv(only_sir_columns = False, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.administrable_iv(only_sir_columns = False, *args, **kwargs))
def not_intrinsic_resistant(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.not_intrinsic_resistant(*args, **kwargs))
def as_ab(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_ab(x, *args, **kwargs))
def is_ab(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_ab(x))
def as_av(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_av(x, *args, **kwargs))
def is_av(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_av(x))
def as_disk(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_disk(x, *args, **kwargs))
def is_disk(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_disk(x))
def as_mic(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_mic(x, *args, **kwargs))
def is_mic(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_mic(x))
def rescale_mic(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.rescale_mic(x, *args, **kwargs))
def as_mo(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_mo(*args, **kwargs))
def is_mo(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_mo(x))
def mo_uncertainties(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_uncertainties(*args, **kwargs))
def mo_renamed(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_renamed(*args, **kwargs))
def mo_failures(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_failures(*args, **kwargs))
def mo_reset_session(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_reset_session(*args, **kwargs))
def mo_cleaning_regex(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_cleaning_regex(*args, **kwargs))
def as_sir(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.as_sir(x, *args, **kwargs))
def is_sir(x):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_sir(x))
def is_sir_eligible(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_sir_eligible(x, *args, **kwargs))
def sir_interpretation_history(clean):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.sir_interpretation_history(clean))
def atc_online_property(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.atc_online_property(*args, **kwargs))
def atc_online_groups(atc_code, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.atc_online_groups(atc_code, *args, **kwargs))
def atc_online_ddd(atc_code, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.atc_online_ddd(atc_code, *args, **kwargs))
def atc_online_ddd_units(atc_code, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.atc_online_ddd_units(atc_code, *args, **kwargs))
def av_from_text(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_from_text(*args, **kwargs))
def av_name(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_name(x, *args, **kwargs))
def av_cid(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_cid(x, *args, **kwargs))
def av_synonyms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_synonyms(x, *args, **kwargs))
def av_tradenames(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_tradenames(x, *args, **kwargs))
def av_group(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_group(x, *args, **kwargs))
def av_atc(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_atc(x, *args, **kwargs))
def av_loinc(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_loinc(x, *args, **kwargs))
def av_ddd(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_ddd(x, *args, **kwargs))
def av_ddd_units(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_ddd_units(x, *args, **kwargs))
def av_info(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_info(x, *args, **kwargs))
def av_url(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_url(x, *args, **kwargs))
def av_property(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.av_property(x, *args, **kwargs))
def availability(tbl, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.availability(tbl, *args, **kwargs))
def bug_drug_combinations(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.bug_drug_combinations(x, *args, **kwargs))
def count_resistant(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_resistant(*args, **kwargs))
def count_susceptible(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_susceptible(*args, **kwargs))
def count_S(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_S(*args, **kwargs))
def count_SI(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_SI(*args, **kwargs))
def count_I(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_I(*args, **kwargs))
def count_IR(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_IR(*args, **kwargs))
def count_R(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_R(*args, **kwargs))
def count_all(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_all(*args, **kwargs))
def n_sir(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.n_sir(*args, **kwargs))
def count_df(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.count_df(*args, **kwargs))
def custom_eucast_rules(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.custom_eucast_rules(*args, **kwargs))
def eucast_rules(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.eucast_rules(*args, **kwargs))
def eucast_dosage(ab, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.eucast_dosage(ab, *args, **kwargs))
def export_ncbi_biosample(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.export_ncbi_biosample(*args, **kwargs))
def first_isolate(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.first_isolate(*args, **kwargs))
def filter_first_isolate(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.filter_first_isolate(*args, **kwargs))
def g_test(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.g_test(x, *args, **kwargs))
def is_new_episode(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.is_new_episode(x, *args, **kwargs))
def ggplot_pca(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ggplot_pca(*args, **kwargs))
def ggplot_sir(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ggplot_sir(*args, **kwargs))
def geom_sir(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.geom_sir(*args, **kwargs))
def theme_sir(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.theme_sir(*args, **kwargs))
def labels_sir_count(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.labels_sir_count(*args, **kwargs))
def guess_ab_col(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.guess_ab_col(*args, **kwargs))
def italicise_taxonomy(string, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.italicise_taxonomy(string, *args, **kwargs))
def italicize_taxonomy(string, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.italicize_taxonomy(string, *args, **kwargs))
def inner_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.inner_join_microorganisms(x, *args, **kwargs))
def left_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.left_join_microorganisms(x, *args, **kwargs))
def right_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.right_join_microorganisms(x, *args, **kwargs))
def full_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.full_join_microorganisms(x, *args, **kwargs))
def semi_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.semi_join_microorganisms(x, *args, **kwargs))
def anti_join_microorganisms(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.anti_join_microorganisms(x, *args, **kwargs))
def key_antimicrobials(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.key_antimicrobials(*args, **kwargs))
def all_antimicrobials(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.all_antimicrobials(x = None, *args, **kwargs))
def antimicrobials_equal(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.antimicrobials_equal(*args, **kwargs))
def kurtosis(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.kurtosis(x, *args, **kwargs))
def like(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.like(x, *args, **kwargs))
def mdro(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mdro(*args, **kwargs))
def custom_mdro_guideline(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.custom_mdro_guideline(*args, **kwargs))
def brmo(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.brmo(x = None, *args, **kwargs))
def mrgn(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mrgn(x = None, *args, **kwargs))
def mdr_tb(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mdr_tb(x = None, *args, **kwargs))
def mdr_cmi2012(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mdr_cmi2012(x = None, *args, **kwargs))
def eucast_exceptional_phenotypes(x = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.eucast_exceptional_phenotypes(x = None, *args, **kwargs))
def mean_amr_distance(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mean_amr_distance(x, *args, **kwargs))
def amr_distance_from_row(amr_distance, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.amr_distance_from_row(amr_distance, *args, **kwargs))
def mo_matching_score(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_matching_score(x, *args, **kwargs))
def mo_name(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_name(*args, **kwargs))
def mo_fullname(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_fullname(*args, **kwargs))
def mo_shortname(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_shortname(*args, **kwargs))
def mo_subspecies(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_subspecies(*args, **kwargs))
def mo_species(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_species(*args, **kwargs))
def mo_genus(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_genus(*args, **kwargs))
def mo_family(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_family(*args, **kwargs))
def mo_order(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_order(*args, **kwargs))
def mo_class(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_class(*args, **kwargs))
def mo_phylum(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_phylum(*args, **kwargs))
def mo_kingdom(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_kingdom(*args, **kwargs))
def mo_domain(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_domain(*args, **kwargs))
def mo_type(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_type(*args, **kwargs))
def mo_status(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_status(*args, **kwargs))
def mo_pathogenicity(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_pathogenicity(*args, **kwargs))
def mo_gramstain(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_gramstain(*args, **kwargs))
def mo_is_gram_negative(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_is_gram_negative(*args, **kwargs))
def mo_is_gram_positive(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_is_gram_positive(*args, **kwargs))
def mo_is_yeast(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_is_yeast(*args, **kwargs))
def mo_is_intrinsic_resistant(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_is_intrinsic_resistant(*args, **kwargs))
def mo_oxygen_tolerance(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_oxygen_tolerance(*args, **kwargs))
def mo_is_anaerobic(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_is_anaerobic(*args, **kwargs))
def mo_snomed(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_snomed(*args, **kwargs))
def mo_ref(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_ref(*args, **kwargs))
def mo_authors(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_authors(*args, **kwargs))
def mo_year(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_year(*args, **kwargs))
def mo_lpsn(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_lpsn(*args, **kwargs))
def mo_mycobank(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_mycobank(*args, **kwargs))
def mo_gbif(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_gbif(*args, **kwargs))
def mo_rank(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_rank(*args, **kwargs))
def mo_taxonomy(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_taxonomy(*args, **kwargs))
def mo_synonyms(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_synonyms(*args, **kwargs))
def mo_current(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_current(x, *args, **kwargs))
def mo_group_members(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_group_members(*args, **kwargs))
def mo_info(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_info(*args, **kwargs))
def mo_url(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_url(*args, **kwargs))
def mo_property(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.mo_property(*args, **kwargs))
def pca(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.pca(*args, **kwargs))
def resistance(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.resistance(*args, **kwargs))
def susceptibility(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.susceptibility(*args, **kwargs))
def sir_confidence_interval(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.sir_confidence_interval(*args, **kwargs))
def proportion_R(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_R(*args, **kwargs))
def proportion_IR(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_IR(*args, **kwargs))
def proportion_I(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_I(*args, **kwargs))
def proportion_SI(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_SI(*args, **kwargs))
def proportion_S(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_S(*args, **kwargs))
def proportion_df(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.proportion_df(*args, **kwargs))
def sir_df(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.sir_df(*args, **kwargs))
def random_mic(size = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.random_mic(size = None, *args, **kwargs))
def random_disk(size = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.random_disk(size = None, *args, **kwargs))
def random_sir(size = None, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.random_sir(size = None, *args, **kwargs))
def resistance_predict(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.resistance_predict(*args, **kwargs))
def sir_predict(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.sir_predict(*args, **kwargs))
def ggplot_sir_predict(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.ggplot_sir_predict(*args, **kwargs))
def skewness(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.skewness(x, *args, **kwargs))
def reset_AMR_locale(*args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.reset_AMR_locale(*args, **kwargs))
def translate_AMR(x, *args, **kwargs):
    """See our website of the R package for the manual: https://msberends.github.io/AMR/index.html"""
    return convert_to_python(amr_r.translate_AMR(x, *args, **kwargs))

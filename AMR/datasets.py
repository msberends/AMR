import pandas as pd
from rpy2 import robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter, numpy2ri, pandas2ri

from ._engine import ensure_amr, restore_sink

_cache = {}
_loaded_source = None

def _load_datasets(source="cran"):
    """Load all AMR datasets into the module cache."""
    global _loaded_source
    
    if _cache and _loaded_source == source:
        return
    
    if _cache and _loaded_source != source:
        _cache.clear()
    
    ensure_amr(source)
    
    with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
        _cache['example_isolates'] = _load_example_isolates()
        _cache['microorganisms'] = robjects.r(
            'AMR::microorganisms[, !sapply(AMR::microorganisms, is.list)]')
        _cache['antimicrobials'] = robjects.r(
            'AMR::antimicrobials[, !sapply(AMR::antimicrobials, is.list)]')
        _cache['clinical_breakpoints'] = robjects.r(
            'AMR::clinical_breakpoints[, !sapply(AMR::clinical_breakpoints, is.list)]')
    
    restore_sink()
    _loaded_source = source

def _load_example_isolates():
    df = robjects.r('''
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
    df['date'] = pd.to_datetime(df['date'])
    return df

def get(name, source="cran"):
    """Retrieve a dataset by name, installing AMR if needed."""
    _load_datasets(source)
    return _cache[name]

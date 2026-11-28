import sys

_DATASETS = frozenset({
    'example_isolates', 'microorganisms',
    'antimicrobials', 'clinical_breakpoints'
})

class _BetaModule(type(sys.modules[__name__])):
    """Lazy-loading module: installs AMR from GitHub on first access."""
    
    def __getattr__(self, name):
        if name in _DATASETS:
            from .datasets import get
            return get(name, source="github")
        try:
            from . import functions
            return getattr(functions, name)
        except AttributeError:
            raise AttributeError(
                f"module 'AMR.beta' has no attribute '{name}'")

sys.modules[__name__].__class__ = _BetaModule

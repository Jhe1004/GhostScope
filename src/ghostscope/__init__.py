"""GhostScope: conditional evidence with an explicit candidate universe."""
__version__ = "1.0.0"

def audit_trees(*args, **kwargs):
    from .analysis import audit_trees as run
    return run(*args, **kwargs)

def audit_chains(*args, **kwargs):
    from .mcmc_engine import audit_chains as run
    return run(*args, **kwargs)

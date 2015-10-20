"""
Package containing shared third-party code modules included here to reduce the need for external
dependencies when only small parts of external code are used.
"""
def midas_available():
    try:
        import MIDAS
        return True
    except ImportError:
        return False

__all__ = []

if midas_available():
    __all__ = ["omsi_score_compounds"]
    from omsi.analysis.compound_stats.omsi_score_compounds import omsi_score_compounds

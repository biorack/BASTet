"""
Package for studies of chemical compounds
"""
def midas_available():
    try:
        import MIDAS
        return True
    except ImportError:
        return False

__all__ = []

if midas_available():
    __all__ = ["omsi_midas"]
    from omsi.analysis.compound_stats.omsi_midas import omsi_midas

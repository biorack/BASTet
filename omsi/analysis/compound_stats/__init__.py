"""
Package containing shared third-party code modules included here to reduce the need for external
dependencies when only small parts of external code are used.
"""
def midas_available():
    """
    Function used to check whether MIDAS can be loaded.

    :return: Boolean indicating whether MIDAS is available
    """
    try:
        import MIDAS
        return True
    except ImportError:
        return False

def pactolus_available():
    """
    Function used to check whether Pactolus can be loaded

    :return: Boolean indicating whether Pactolus is available
    """
    try:
        import pactolus
        return True
    except ImportError:
        return False


__all__ = []

if midas_available():
    __all__ = ["omsi_score_midas"]
    from omsi.analysis.compound_stats.omsi_score_midas import omsi_score_midas

if pactolus_available():
    __all__ = ["omsi_score_pactolus"]
    from omsi.analysis.compound_stats.omsi_score_pactolus import omsi_score_pactolus

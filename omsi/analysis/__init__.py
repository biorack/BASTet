"""
Package containing the base classes that facilitate the integration of new analysis with the BASTet software
stack (e.g, the file format) and collection of specific analysis functionality.
"""

from omsi.analysis.findpeaks import *
from omsi.analysis.multivariate_stats import *
from omsi.analysis.msi_filtering import *
from omsi.analysis.base import *
from omsi.analysis.generic import *
from omsi.analysis.compound_stats import *
import omsi.analysis.findpeaks
import omsi.analysis.multivariate_stats
import omsi.analysis.msi_filtering
import omsi.analysis.compound_stats
# from analysis_views import *
__all__ = ["findpeaks",
           "multivariate_stats",
           "compound_stats",
           "msi_filtering",
           "analysis_data",
           "analysis_base",
           "AnalysisReadyError",
           "analysis_generic"] + findpeaks.__all__ + multivariate_stats.__all__ + msi_filtering.__all__ + compound_stats.__all__


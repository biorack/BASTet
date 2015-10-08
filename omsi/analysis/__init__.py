"""
OpenMSI toolkit analysis package defining the main analysis API as well as
specific MSI analysis algorithms.
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
           "analysis_generic"] + findpeaks.all__ + multivariate_stats.all__ + msi_filtering.all__


"""
Multivariate statistics analysis
"""


def kmeans_available():
    try:
        from scipy.cluster.vq import kmeans, vq
        return True
    except ImportError:
        return False

__all__ = ["omsi_nmf", "omsi_cx"]
from omsi_nmf import *
from omsi_cx import *
if kmeans_available():
    __all__ += ["omsi_kmeans"]
    from omsi_kmeans import *

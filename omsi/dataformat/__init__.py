"""
Package for implementation and specification of file formats.
"""

def mzml_available():
    """
    Check whether additional optional libraries are
    available that are required for the use of the
    mzml file reader.
    """
    try:
        import pyteomics
        import lxml
        return True
    except ImportError:
        return False

def imzml_available():
    """
    Check whether addtional optional libraries are
    available that are required for the use of teh
    imzml file reader.
    """
    try:
        from pyimzml.ImzMLParser import ImzMLParser
        return True
    except ImportError:
        return False

__all__ = ["img_file", "bruckerflex_file", "omsi_file", "file_reader_base"]
if mzml_available():
    __all__.append("mzml_file")
    __all__.append("xmassmzml_file")
if imzml_available():
    __all__.append("imzml_file")


# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A package for computing orbits and astronomical observables for binary stars,
exoplanets, and other gravitational two-body systems.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__minimum_python_version__ = "3.6"

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("packagename does not support Python < {}"
                                 .format(__minimum_python_version__))

if not _ASTROPY_SETUP_:
    from .anomaly import *
    from .barycenter import *
    from .bary_trends import *
    from .elements import *
    from .orbit import *
    from .reference_plane import *
    from .transforms import *
    from .units import *

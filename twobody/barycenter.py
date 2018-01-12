# Third-party
from astropy.constants import G
from astropy.time import Time
import astropy.units as u
import numpy as np

# Project
from .utils import _parse_time

__all__ = ['Barycenter']


class Barycenter:

    def __init__(self, origin, t0=None):
        """The location of the barycenter at the specified epoch.

        Parameters
        ----------
        origin : `~astropy.coordinates.SkyCoord` or frame instance
            TODO
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch for the location of the barycenter. If a number is
            passed in, it is assumed to be a solar system barycentric modified
            julian date (BMJD). The default is J2000 if not specified.
        """
        self.origin = origin

        if t0 is None:
            t0 = Time('J2000')
        self.t0 = _parse_time(t0)

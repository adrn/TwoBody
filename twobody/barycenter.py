# Third-party
from astropy.constants import G
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy import pi

# Project
from .utils import _parse_time


__all__ = ['Barycenter']


class Barycenter:

    def __init__(self, coords, t0):
        self.coords = coords
        self.t0 = _parse_time(t0)


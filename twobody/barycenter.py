# Standard library
import warnings

# Third-party
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import numpy as np

# Project
from .elements import _parse_time

__all__ = ['Barycenter']


class Barycenter:

    def __init__(self, origin=None, t0=None, radial_velocity=None):
        """The location of the barycenter at the specified epoch.

        This class is used to specify a reference location and velocity of the
        barycenter at a given epoch. This is used for computing astrometric
        orbits, and for specifying the barycentric motion.

        In the future, this class will contain more functionality that supports,
        e.g., acceleration of the barycenter from a third body.

        Parameters
        ----------
        origin : `~astropy.coordinates.SkyCoord` or frame instance
            The sky position of and distance to the barycenter at the specified
            epoch, ``t0``. If only a barycentric radial velocity is known or
            needed, pass in ``radial_velocity``.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch for the location of the barycenter. If a number is
            passed in, it is assumed to be a solar system barycentric modified
            julian date (BMJD). The default is J2000 if not specified.
        radial_velocity : quantity_like (optional)
            If the coordinates or distance to the barycenter are not known or
            needed, you can just pass the radial velocity of the barycenter
            instead.

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.coordinates import SkyCoord
        >>> b = Barycenter(radial_velocity=100*u.km/u.s)
        >>> b = Barycenter(origin=SkyCoord(ra=150*u.deg, dec=37*u.deg,
        ...                                distance=100*u.pc)) # assumes J2000
        >>> from astropy.time import Time
        >>> b = Barycenter(origin=SkyCoord(ra=150*u.deg, dec=37*u.deg,
        ...                                distance=100*u.pc),
        ...                t0=Time(58912.481293, format='mjd'))

        """

        if origin is not None and radial_velocity is not None:
            raise ValueError("Pass one of: origin or radial_velocity, not both")

        if origin is None and radial_velocity is None:
            raise ValueError("You must pass one of origin or radial_velocity")

        if origin is None:
            origin = coord.ICRS(ra=0*u.deg, dec=0*u.deg, distance=np.nan*u.pc,
                                radial_velocity=radial_velocity)

        if (not isinstance(origin, coord.SkyCoord) and
                not isinstance(origin, coord.BaseCoordinateFrame)):
            raise TypeError("origin must be a SkyCoord or BaseCoordinateFrame "
                            "subclass instance, i.e. a valid "
                            "astropy.coordinates object.")

        if not origin.isscalar:
            raise ValueError("origin must be a scalar coordinate, not "
                             "array-valued: {0}".format(origin))

        if origin.distance.unit == u.dimensionless_unscaled:
            warnings.warn("Barycenter origin has no associated distance. Are "
                          "you sure you know what you're doing?",
                          UserWarning)

        self.origin = origin

        if t0 is None:
            t0 = Time('J2000')
        self.t0 = _parse_time(t0)

    def __repr__(self):
        fmt_str = "<Barycenter: origin={0.origin}, epoch={0.t0}>".format(self)
        return fmt_str

    def __str__(self):
        return repr(self)

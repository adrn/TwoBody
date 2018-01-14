# Standard library
import warnings

# Third-party
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u

# Project
from .elements import _parse_time

__all__ = ['Barycenter']


class Barycenter:

    def __init__(self, origin, t0=None):
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
            epoch, ``t0``.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch for the location of the barycenter. If a number is
            passed in, it is assumed to be a solar system barycentric modified
            julian date (BMJD). The default is J2000 if not specified.

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.coordinates import SkyCoord
        >>> b = Barycenter(origin=SkyCoord(ra=150*u.deg, dec=37*u.deg,
        ...                                distance=100*u.pc)) # assumes J2000
        >>> from astropy.time import Time
        >>> b = Barycenter(origin=SkyCoord(ra=150*u.deg, dec=37*u.deg,
        ...                                distance=100*u.pc),
        ...                t0=Time(58912.481293, format='mjd'))

        """

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

    def __str__(self):
        return repr(self)

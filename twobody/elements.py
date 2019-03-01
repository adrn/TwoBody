# Standard library
import abc
import warnings

# Third-party
import astropy.coordinates as coord
from astropy.constants import G, c
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy import pi

# Project
from .transforms import a_m_to_P, P_m_to_a, PeKi_to_a
from .units import UnitSystem

__all__ = ['OrbitalElements', 'KeplerElements', 'TwoBodyKeplerElements']


def _parse_time(t):
    if not isinstance(t, Time):
        return Time(t, format='mjd', scale='tcb')
    return t


class ElementsMeta(abc.ABCMeta):

    def __new__(mcls, name, bases, members):

        if 'names' not in members:
            raise ValueError('OrbitalElements subclasses must contain a '
                             'defined class attribute "names" that specified '
                             'the string names of the elements.')

        for name_ in members['names']:
            mcls.readonly_prop_factory(members, name_)

        return super().__new__(mcls, name, bases, members)

    @staticmethod
    def readonly_prop_factory(members, attr_name):
        def getter(self):
            return self.units.decompose(getattr(self, '_' + attr_name))
        members[attr_name] = property(getter)


class OrbitalElements(metaclass=ElementsMeta):
    """
    Subclasses must define the class attribute ``.default_units`` to be a
    ``UnitSystem`` instance.
    """

    names = []

    def __init__(self, units):

        # Make sure the units specified are a UnitSystem instance
        if units is None:
            units = self.default_units

        elif units is not None and not isinstance(units, UnitSystem):
            units = UnitSystem(*units)

        self.units = units

        # Now make sure all element name attributes have been set:
        for name in self.names:
            if not hasattr(self, '_' + name):
                raise AttributeError('Invalid class definition!')


class BaseKeplerElements(OrbitalElements):
    default_units = UnitSystem(u.au, u.day, u.Msun, u.degree, u.km / u.s)
    names = ['P', 'a', 'e', 'omega', 'i', 'Omega', 'M0']

    def __init__(self, P=None, a=None,
                 e=0, omega=None, i=None, Omega=None,
                 M0=None, t0=None, units=None):
        """Base class for Keplerian orbital elements.

        Parameters
        ----------
        P : quantity_like [time]
            Orbital period.
        a : quantity_like [length] (optional)
            Semi-major axis. If unspecified, computed orbits will be unscaled.
        e : numeric (optional)
            Orbital eccentricity. Default is circular, ``e=0``.
        omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Argument of pericenter.
        i : quantity_like, `~astropy.coordinates.Angle` [angle]
            Inclination of the orbit.
        Omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Longitude of the ascending node.
        M0 : quantity_like, `~astropy.coordinates.Angle` [angle] (optional)
            Mean anomaly at epoch ``t0``. Default is 0º if not specified.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch. If a number is passed in, it is assumed to be
            a solar system barycentric modified julian date (BMJD). The default
            is J2000 if not specified.
        units : `~twobody.units.UnitSystem`, iterable (optional)
            The unit system to represent quantities in. The default unit system
            is accessible as `KeplerElements.default_units`.

        """

        if M0 is None:
            # Default phase at reference epoch is 0º
            M0 = 0 * u.degree

        if t0 is None:
            # Default reference epoch is J2000
            t0 = Time('J2000')

        t0 = _parse_time(t0)

        # Now check that required elements are defined:
        _required = ['P', 'omega']
        for name in _required:
            if eval(name) is None:
                raise ValueError("You must specify {0}.".format(name))

        # Value validation:
        if P < 0 * u.day:
            raise ValueError("Period `P` must be positive.")

        if a is not None and a < 0 * u.au:
            raise ValueError("Semi-major axis `a` must be positive.")

        if e < 0 or e >= 1:
            raise ValueError("Eccentricity `e` must be: 0 <= e < 1")

        if i is not None and (i < 0*u.deg or i > 180 * u.deg):
            raise ValueError("Inclination `i` must be between 0º and 180º, you "
                             "passed in i={:.3f}".format(i.to(u.degree)))

        if i is None:
            i = np.nan * u.deg

        if Omega is None:
            Omega = np.nan * u.deg

        # Set object attributes, but make them read-only
        self._a = a if a is not None else 1. * u.dimensionless_unscaled
        self._P = P
        self._e = float(e) * u.dimensionless_unscaled
        self._omega = coord.Angle(omega).wrap_at(360 * u.deg)
        self._i = coord.Angle(i)
        self._Omega = coord.Angle(Omega).wrap_at(360 * u.deg)
        self._M0 = coord.Angle(M0)
        self.t0 = t0

        # Must happen at the end because it validates that all element names
        # have been set properly:
        super().__init__(units=units)


class KeplerElements(BaseKeplerElements):
    names = ['P', 'a', 'e', 'omega', 'i', 'Omega', 'M0']

    @u.quantity_input(P=u.year, a=u.au, K=u.km/u.s,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, P=None, a=None, K=None,
                 e=0, omega=None, i=None, Omega=None,
                 M0=None, t0=None, units=None):
        """Keplerian orbital elements for a single orbit.

        The elements are assumed to be relative to an inertial frame, typically
        the barycenter of a two-body system.

        Parameters
        ----------
        P : quantity_like [time]
            Orbital period.
        a : quantity_like [length] (optional)
            Semi-major axis. Specify this OR the semi-amplitude ``K``, but not
            both. If unspecified, computed orbits will be unscaled.
        K : quantity_like [speed] (optional)
            Velocity semi-amplitudes. Specify this OR the semi-major axis ``a``,
            but not both. If unspecified, computed orbits will be unscaled.
        e : numeric (optional)
            Orbital eccentricity. Default is circular, ``e=0``.
        omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Argument of pericenter.
        i : quantity_like, `~astropy.coordinates.Angle` [angle]
            Inclination of the orbit.
        Omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Longitude of the ascending node.
        M0 : quantity_like, `~astropy.coordinates.Angle` [angle] (optional)
            Mean anomaly at epoch ``t0``. Default is 0º if not specified.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch. If a number is passed in, it is assumed to be
            a solar system barycentric modified julian date (BMJD). The default
            is J2000 if not specified.
        units : `~twobody.units.UnitSystem`, iterable (optional)
            The unit system to represent quantities in. The default unit system
            is accessible as `KeplerElements.default_units`.

        """
        if K is not None and a is not None:
            raise ValueError("Both `K` and `a` were specified, but can only "
                             "accept one or the other.")

        super().__init__(P=P, a=a, e=e, omega=omega, i=i, Omega=Omega,
                         M0=M0, t0=t0, units=units)

        # TODO: this is a bit messy because we double-initialize! The first
        # time to validate, the second to use the proper "a"
        if K is not None:
            a = PeKi_to_a(P, e, K, i)
            super().__init__(P=P, a=a, e=e, omega=omega, i=i, Omega=Omega,
                             M0=M0, t0=t0, units=units)

        if self.K.unit.physical_type == 'speed':
            # Now we do a quick validation to make sure we're not relativistic...
            _K_c = (self.K / c).decompose()
            if _K_c > 1E-2:
                warnings.warn('Velocity semiamplitude is large enough that '
                              'relativistic effects are important: K/c = {:.2f}'
                              ' but are not currently accounted for in TwoBody'
                              .format(_K_c), RuntimeWarning)

    @property
    def K(self):
        """Velocity semi-amplitude."""
        K = 2*pi * self.a * np.sin(self.i) / (self.P * np.sqrt(1 - self.e**2))
        return self.units.decompose(K)

    @property
    def m_f(self):
        """Binary mass function."""
        mf_circ = self.P * self.K**3 / (2*pi * G)
        return self.units.decompose(mf_circ) * (1 - self.e**2)**1.5

    # Python builtins
    def __repr__(self):
        return ("<KeplerElements [P={:.2f}, a={:.2f}, e={:.2f}, "
                "ω={:.2f}, i={:.2f}, Ω={:.2f}]>"
                .format(self.P, self.a, self.e, self.omega, self.i, self.Omega))


class TwoBodyKeplerElements(BaseKeplerElements):
    names = ['P', 'a', 'e', 'm1', 'm2', 'omega', 'i', 'Omega', 'M0']

    @u.quantity_input(P=u.year, a=u.au, m1=u.Msun, m2=u.Msun,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, P=None, a=None, m1=None, m2=None,
                 e=0, omega=None, i=None, Omega=None,
                 M0=None, t0=None, units=None):
        """Keplerian orbital elements for a two-body system.

        You can either specify period, ``P`` and the masses ``m1``, ``m2``, or
        the separation ``a`` and the two masses.

        Parameters
        ----------
        P : quantity_like [time]
            Orbital period.
        a : quantity_like [length]
            Semi-major axis of the fictitious particle
        e : numeric (optional)
            Orbital eccentricity. Default is circular, ``e=0``.
        omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Argument of pericenter.
        i : quantity_like, `~astropy.coordinates.Angle` [angle]
            Inclination of the orbit.
        Omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Longitude of the ascending node of the fictitious particle.
        M0 : quantity_like, `~astropy.coordinates.Angle` [angle] (optional)
            Mean anomaly at epoch ``t0``. Default is 0º if not specified.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch. Default is J2000 if not specified.
        units : `~twobody.units.UnitSystem`, iterable (optional)
            The unit system to represent quantities in. The default unit system
            is accessible as `KeplerElements.default_units`.

        """

        if m1 is None or m2 is None:
            raise ValueError("You must specify m1 and m2.")

        if P is not None and a is not None:
            raise ValueError("You can only specify one of period `P` or "
                             "semi-major axis `a`.")

        if P is None:
            P = a_m_to_P(a, m1 + m2)

        if a is None:
            a = P_m_to_a(P, m1 + m2)

        self._m1 = m1
        self._m2 = m2
        self.m_tot = m1 + m2

        # values are validated
        super().__init__(a=a, P=P,
                         e=e, omega=omega, i=i, Omega=omega,
                         M0=M0, t0=t0, units=units)

    def get_body(self, num):
        """Get the orbital elements for the specified body.

        Parameters
        ----------
        num : str ('1' or '2')
            Get the orbital elements of the primary `'1'` or secondary `'2'`.

        Returns
        -------
        elements : `twobody.KeplerElements`

        """
        num = str(num)

        if num == '1':
            a = self.m2 / self.m_tot * self.a
            omega = self.omega

        elif num == '2':
            a = self.m1 / self.m_tot * self.a
            omega = self.omega + np.pi*u.radian

        else:
            raise ValueError("Invalid input '{0}' - must be '1' or '2'"
                             .format(num))

        return KeplerElements(a=a, P=self.P,
                              e=self.e, omega=omega, i=self.i, Omega=self.Omega,
                              M0=self.M0, t0=self.t0)

    @property
    def primary(self):
        """Shorthand for ``.get_body('1')``. Returns the orbital elements for
        the primary body, which corresponds to ``m1`` (not the more massive of
        the two).
        """
        return self.get_body('1')

    @property
    def secondary(self):
        """Shorthand for ``.get_body('2')``. Returns the orbital elements for
        the secondary body, which corresponds to ``m2`` (not the less massive of
        the two).
        """
        return self.get_body('2')

    # Python builtins
    def __repr__(self):
        return ("<TwoBodyKeplerElements [m1={:.2f}, m2={:.2f}, "
                "P={:.2f}, a={:.2f}, e={:.2f}, "
                "ω={:.2f}, i={:.2f}, Ω={:.2f}]>"
                .format(self.m1, self.m2, self.P, self.a, self.e,
                        self.omega, self.i, self.Omega))

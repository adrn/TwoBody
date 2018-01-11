# Third-party
from astropy.constants import G
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy import pi

# Project
from .transforms import a_m_to_P, P_m_to_a
from .units import UnitSystem

__all__ = ['OrbitalElements', 'KeplerElements', 'TwoBodyKeplerElements']


class OrbitalElements:

    def __init__(self, units):

        # make sure the units specified are a UnitSystem instance
        if units is not None and not isinstance(units, UnitSystem):
            units = UnitSystem(*units)

        self.units = units

class KeplerElements(OrbitalElements):

    default_units = UnitSystem(u.au, u.day, u.Msun, u.degree, u.km/u.s)

    @u.quantity_input(a=u.au, P=u.year,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, a=None, P=None,
                 e=0, omega=None, i=None, Omega=None,
                 M0=None, t0=None, units=None):
        """Class for representing Keplerian orbital elements.

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
            Reference epoch. Default is J2000 if not specified.
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

        if not isinstance(t0, Time):
            # If a number is specified, assume it is Barycentric MJD
            t0 = Time(t0, format='mjd', scale='tcb')

        # Default unit system:
        if units is None:
            units = self.default_units
        super().__init__(units=units)

        # Now check that required elements are defined:
        _required = ['P', 'omega', 'i', 'Omega']
        for name in _required:
            if eval(name) is None:
                raise ValueError("You must specify {0}.".format(name))

        if i < 0*u.deg or i > 180*u.deg:
            raise ValueError("Inclination must be between 0º and 180º, you "
                             "passed in i={:.3f}".format(i.to(u.degree)))

        # Set object attributes
        self.a = a if a is not None else 1.
        self.P = P
        self.e = float(e)
        self.omega = omega
        self.i = i
        self.Omega = Omega
        self.M0 = M0
        self.t0 = t0

    @property
    def K(self):
        """Velocity semi-amplitude."""
        return 2*pi * self.a * np.sin(self.i) / (self.P * np.sqrt(1-self.e**2))

    @property
    def m_f(self):
        """Binary mass function."""
        return self.P * self.K**3 / (2*pi * G)

    # Python builtins
    def __repr__(self):
        return ("<KeplerElements [P={:.2f}, a={:.2f}, e={:.2f}, "
                "ω={:.2f}, i={:.2f}, Ω={:.2f}]>"
                .format(self.P, self.a, self.e, self.omega, self.i, self.Omega))


# TODO: be very explicit. Are we specifying the elements of one of the bodies,
# or of the fictitious body? Or allow user to specify?
class TwoBodyKeplerElements(KeplerElements):

    @u.quantity_input(a=u.au, P=u.year, m1=u.Msun, m2=u.Msun,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, a=None, P=None, m1=None, m2=None,
                 e=None, omega=None, i=None, Omega=None,
                 M0=None, t0=None):

        if m1 is None or m2 is None:
            raise ValueError("You must specify m1 and m2.")

        if P is None:
            P = a_m_to_P(a, m1+m2).to(u.day)

        if a is None:
            a = P_m_to_a(P, m1+m2).to(u.au)

        super().__init__(a=a, P=P,
                         e=e, omega=omega, i=i, Omega=omega,
                         M0=M0, t0=t0)

        self.m1 = m1
        self.m2 = m2
        self.m_tot = self.m1 + self.m2

    def get_component(self, num):
        """TODO
        """
        num = str(num)

        if num == '1':
            a = self.m2 / self.m_tot * self.a
            omega = self.omega
            Omega = self.Omega

        elif num == '2':
            a = self.m1 / self.m_tot * self.a
            omega = self.omega + np.pi*u.radian
            Omega = self.Omega + np.pi*u.radian

        else:
            raise ValueError("Invalid input '{0}' - must be '1' or '2'"
                             .format(num))

        return KeplerElements(a=a, P=self.P,
                              e=self.e, omega=omega, i=self.i, Omega=Omega,
                              M0=self.M0, t0=self.t0)

    @property
    def primary(self):
        return self.get_component('1')

    @property
    def secondary(self):
        return self.get_component('2')

    # Python builtins
    def __repr__(self):
        return ("<TwoBodyKeplerElements [m1={:.2f}, m2={:.2f}, "
                "P={:.2f}, a={:.2f}, e={:.2f}, "
                "ω={:.2f}, i={:.2f}, Ω={:.2f}]>"
                .format(self.m1, self.m2, self.P, self.a, self.e,
                        self.omega, self.i, self.Omega))

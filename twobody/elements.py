# Third-party
from astropy.constants import G
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy import pi

# Project
from .anomaly import (eccentric_anomaly_from_mean_anomaly,
                      true_anomaly_from_eccentric_anomaly)
from .wrap import cy_rv_from_elements
from .utils import ArrayProcessor
from .transforms import a_P_to_m, a_m_to_P, P_m_to_a


__all__ = ['OrbitalElements', 'KeplerElements', 'TwoBodyKeplerElements']


class OrbitalElements:

    def __init__(self):
        pass


class KeplerElements(OrbitalElements):

    @u.quantity_input(a=u.au, P=u.year,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, a=None, P=None,
                 e=None, omega=None, i=None, Omega=None,
                 M0=None, t0=None):
        """Class for representing Keplerian orbital elements.

        Parameters
        ----------
        a : quantity_like [length]
            Semi-major axis.
        P : quantity_like [time]
            Orbital period.

        e : numeric
            Orbital eccentricity.
        omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Argument of pericenter.
        i : quantity_like, `~astropy.coordinates.Angle` [angle]
            Inclination of the orbit.
        Omega : quantity_like, `~astropy.coordinates.Angle` [angle]
            Longitude of the ascending node.

        M0 : quantity_like, `~astropy.coordinates.Angle` [angle] (optional)
            Mean anomaly at epoch ``t0``. Assumed to be 0º if not specified.
        t0 : numeric, `~astropy.coordinates.Time` (optional)
            Reference epoch. Assumed to be J2000 if not specified.
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

        # Now check that required things exist:
        _required = ['P', 'e', 'omega', 'i', 'Omega']
        for name in _required:
            if eval(name) is None:
                raise ValueError("You must specify {0} or enough elements to "
                                 "compute {0}.".format(name))

        # TODO: wrap angles

        # Set object attributes
        self.a = a
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


class TwoBodyKeplerElements(KeplerElements):

    @u.quantity_input(a=u.au, P=u.year, m1=u.Msun, m2=u.Msun,
                      omega=u.deg, i=u.deg, Omega=u.deg, M0=u.deg)
    def __init__(self, *, a=None, P=None, m1=None, m2=None,
                 e=None, omega=None, i=None, Omega=None,
                 M0=None, t0=None):

        if m1 is None or m2 is None:
            raise ValueError("You must specify m1 and m2.")

        if P is None:
            P = a_m_to_P(a, m1+m2)

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
        return self.get_component('1')

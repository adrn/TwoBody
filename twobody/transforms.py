# Third-party
import astropy.units as u
from astropy.utils.misc import isiterable
from astropy.constants import G
import numpy as np

# Project
from .utils import format_doc


__all__ = ['a_P_to_m', 'a_m_to_P', 'P_m_to_a', 'get_m2_min']


doc_a = """a : quantity_like [length]
        Semi-major axis.
"""

doc_P = """P : quantity_like [time]
        Orbital period.
"""

doc_m = """m : quantity_like [mass]
        Total mass.
"""


@u.quantity_input(a=u.au, P=u.day)
@format_doc("""{__doc__}""", a=doc_a, P=doc_P)
def a_P_to_m(a, P):
    """Compute the total mass given the semi-major axis and period.

    Parameters
    ----------
    {a}
    {P}
    """
    return a**3 / G * (2*np.pi/P)**2


@u.quantity_input(a=u.au, m=u.Msun)
@format_doc("""{__doc__}""", a=doc_a, m=doc_m)
def a_m_to_P(a, m):
    """Compute the orbital period given the semi-major axis and total mass.

    Parameters
    ----------
    {a}
    {m}
    """
    return 2*np.pi * np.sqrt(a**3 / (G * m))


@u.quantity_input(P=u.day, m=u.Msun)
@format_doc("""{__doc__}""", P=doc_P, m=doc_m)
def P_m_to_a(P, m):
    """Compute the semi-major axis given the orbital period and total mass.

    Parameters
    ----------
    {P}
    {m}
    """
    return np.cbrt(G * m * (P/(2*np.pi))**2)


@u.quantity_input(P=u.day, K=u.km/u.s, i=u.deg)
@format_doc("""{__doc__}""", P=doc_P)
def PeKi_to_a(P, e, K, i=None):
    """Compute the semi-major axis given the orbital period, eccentricity,
    semi-amplitude, and inclination. If you don't have a measured inclination,
    use the default i=90º and the returned semi-major axis will be a*sin(i)
    instead.

    Parameters
    ----------
    {P}
    e : numeric
        Eccentricity.
    K : quantity_like [speed]
        Velocity semi-amplitude.
    i : quantity_like [angle]
        Inclination.

    Returns
    -------
    a : `~astropy.units.Quantity` [length]
        Semi-major axis.
    """
    if i is None:
        i = 90*u.deg
    a = K / (2*np.pi * np.sin(i)) * (P * np.sqrt(1 - e**2))
    return a


def _m2_func(m2, m1, sini, mf):
    return (m2*sini)**3 / (m1 + m2)**2 - mf


@u.quantity_input(m1=u.Msun, mf=u.Msun)
def get_m2_min(m1, mf):
    """Compute the minimum companion mass given the primary mass and the binary
    mass function.

    Parameters
    ----------
    m1 : quantity_like [mass]
        Primary mass.
    mf : quantity_like [mass]
        Binary mass function.

    Returns
    -------
    m2_min : `~astropy.units.Quantity` [mass]
        The minimum companion mass.
    """
    from scipy.optimize import root

    mf = mf.to(m1.unit)
    if isiterable(m1) and isiterable(mf):
        m2s = []
        for x, y in zip(m1, mf):
            try:
                res = root(_m2_func, x0=10., args=(x.value, 1., y.value))
                if not res.success:
                    raise RuntimeError('Unsuccessful')
                m2s.append(res.x[0])
            except Exception as e:
                m2s.append(np.nan)
        return m2s * m1.unit

    else:
        res = root(_m2_func, x0=10., args=(m1.value, 1., mf.value))

        if res.success:
            return res.x[0] * m1.unit
        else:
            return np.nan * m1.unit

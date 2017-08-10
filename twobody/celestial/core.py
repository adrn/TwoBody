# Third-party
import astropy.units as u
import numpy as np
au_per_day_m_s = (1*u.m/u.s*u.day).to(u.au).value

from .anomaly import (eccentric_anomaly_from_mean_anomaly,
                      true_anomaly_from_eccentric_anomaly)

__all__ = ['Z_from_elements', 'rv_from_elements']

def Z_from_elements(times, P, K, e, omega, time0):
    """
    Z points towards the observer.

    Parameters
    ----------
    times : array_like [day]
        BJD of observations.
    p : numeric [day]
        Period.
    K : numeric [m/s]
        Velocity semi-amplitude.
    e : numeric
        Eccentricity.
    omega : numeric [radian]
        Perihelion argument parameter from Winn.
    time0 : numeric [day]
        Time of "zeroth" pericenter.

    Returns
    -------
    Z : numeric [AU]
        Line-of-sight position.

    TODO
    ----
    - Doesn't include system Z value (Z offset or Z zeropoint) or option to
      specify this quantity. This is the observer-barycenter distance.

    """
    times = np.array(times)

    dMdt = 2. * np.pi / P
    Ms = (times - time0) * dMdt

    Es = eccentric_anomaly_from_mean_anomaly(Ms, e)
    fs = true_anomaly_from_eccentric_anomaly(Es, e)

    a1sini = K/(2*np.pi) * (P * np.sqrt(1-e**2)) * au_per_day_m_s
    rs = a1sini * (1. - e * np.cos(Es))
    # this is equivalent to:
    # rs = asini * (1. - e**2) / (1 + e*np.cos(fs))

    return rs * np.sin(omega + fs)

def rv_from_elements(times, P, K, e, omega, phi0, anomaly_tol=None):
    """
    Parameters
    ----------
    times : array_like [day]
        Usually: Barycentric MJD of observations. But the epoch (t=0)
        is arbitrary and up to the user to keep track of.
    p : numeric [day]
        Period.
    K : numeric [m/s]
        Velocity semi-amplitude.
    e : numeric
        Eccentricity.
    omega : numeric [radian]
        Argument of periastron.
    phi0 : numeric [radian]
        Phase at pericenter relative to t=0.
    anomaly_tol : numeric, optional
        Tolerance passed to
        `~twobody.celestial.eccentric_anomaly_from_mean_anomaly` for solving for
        the eccentric anomaly. See default value in that function.

    Returns
    -------
    rv : numeric [m/s]
        Relative radial velocity - does not include systemtic velocity!
    """
    times = np.array(times)
    Ms = (2 * np.pi * times / P) - phi0

    kw = dict()
    if anomaly_tol is not None:
        kw['tol'] = anomaly_tol
    Es = eccentric_anomaly_from_mean_anomaly(Ms, e, **kw)
    fs = true_anomaly_from_eccentric_anomaly(Es, e)
    vz = K * (np.cos(omega + fs) + e*np.cos(omega))

    return vz

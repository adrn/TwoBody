# Third-party
import astropy.units as u
import numpy as np
au_per_day_m_s = (1*u.m/u.s*u.day).to(u.au).value

from .anomaly import (eccentric_anomaly_from_mean_anomaly,
                      true_anomaly_from_eccentric_anomaly)
from .wrap import cy_rv_from_elements
from .utils import ArrayProcessor

__all__ = ['z_from_elements', 'rv_from_elements']

def z_from_elements(times, P, K, e, omega, time0):
    """Position of the primary along axis connecting the Barycenter and the
    observer.

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
        Relative line-of-sight position; does not include distance to the
        Barycenter of the system!

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

def rv_from_elements(times, P, K, e, omega, phi0,
                     anomaly_tol=None, anomaly_maxiter=None):
    """Compute the dadial velocity of the primary relative to the Barycenter of
    the system at the specified times.

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
    if anomaly_tol is None:
        anomaly_tol = 1E-10

    if anomaly_maxiter is None:
        anomaly_maxiter = 128

    proc = ArrayProcessor(times)
    t, = proc.prepare_arrays()
    rv = cy_rv_from_elements(t, P, K, e, omega, phi0,
                             anomaly_tol, anomaly_maxiter)
    return np.atleast_1d(proc.prepare_result(rv))

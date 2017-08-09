# Third-party
import numpy as np

from .wrap import (cy_mean_anomaly_from_eccentric_anomaly,
                   cy_eccentric_anomaly_from_mean_anomaly_Newton1)

__all__ = ['mean_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_mean_anomaly_Newton1',
           'true_anomaly_from_eccentric_anomaly',
           'd_eccentric_anomaly_d_mean_anomaly',
           'd_true_anomaly_d_eccentric_anomaly']

def mean_anomaly_from_eccentric_anomaly(E, e):
    """
    Parameters
    ----------
    E : numeric, array_like [radian]
        Eccentric anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    M : numeric, array_like [radian]
        Mean anomaly.
    """
    return cy_mean_anomaly_from_eccentric_anomaly(np.ascontiguousarray(E),
                                                  np.ascontiguousarray(e))

# ------------------------------------------------------------------------------

base_docstring = """
Parameters
----------
M : numeric, array_like [radian]
    Mean anomaly.
e : numeric
    Eccentricity.
tol : numeric, optional
    Numerical tolerance used in iteratively solving for eccentric anomaly.
maxiter : int, optional
    Maximum number of iterations when iteratively solving for eccentric
    anomaly.
{0}

Returns
-------
E : numeric [radian]
    Eccentric anomaly.

Issues
------
- Magic numbers ``tol`` and ``maxiter``
"""

def eccentric_anomaly_from_mean_anomaly(M, e, tol=1E-10, maxiter=128,
                                        method='Newton1'):
    func_name = "eccentric_anomaly_from_mean_anomaly_{0}".format(method)
    func = eval(func_name)
    return func(M, e, tol, maxiter)

eccentric_anomaly_from_mean_anomaly.__doc__ = base_docstring.format(
"""method : str, optional
    The method to use for iterative root-finding for the eccentric anomaly.
    Options are: ``'Newton1'``.""")

def eccentric_anomaly_from_mean_anomaly_Newton1(M, e, tol, maxiter):
    return cy_eccentric_anomaly_from_mean_anomaly_Newton1(
        np.ascontiguousarray(M), np.ascontiguousarray(e), tol, maxiter)
eccentric_anomaly_from_mean_anomaly_Newton1.__doc__ = base_docstring.format('')

# ------------------------------------------------------------------------------

def true_anomaly_from_eccentric_anomaly(E, e):
    """
    Parameters
    ----------
    E : numeric, array_like [radian]
        Eccentric anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    fs : numeric [radian]
        True anomaly.
    """
    cEs, sEs = np.cos(E), np.sin(E)
    fs = np.arccos((cEs - e) / (1.0 - e * cEs))
    return fs * np.sign(np.sin(fs)) * np.sign(sEs)

def d_eccentric_anomaly_d_mean_anomaly(E, e):
    """
    Parameters
    ----------
    E : numeric, array_like [radian]
        Eccentric anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    dE_dM : numeric
        Derivatives of eccentric anomaly w.r.t. mean anomaly.
    """
    return 1. / (1. - e * np.cos(E))

def d_true_anomaly_d_eccentric_anomaly(E, f, e):
    """
    Parameters
    ----------
    E : numeric, array_like [radian]
        Eccentric anomaly.
    f : numeric, array_like [radian]
        True anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    df_dE : numeric
        Derivatives of true anomaly w.r.t. eccentric anomaly.

    Issues
    ------
    - Insane assert statement.
    """
    cfs, sfs = np.cos(f), np.sin(f)
    cEs, sEs = np.cos(E), np.sin(E)
    # assert np.allclose(cEs, (e + cfs) / (1. + e * cfs))
    return (sEs / sfs) * (1. + e * cfs) / (1. - e * cEs)

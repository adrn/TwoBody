# Project
from .wrap import (cy_mean_anomaly_from_eccentric_anomaly,
                   cy_true_anomaly_from_eccentric_anomaly,
                   cy_eccentric_anomaly_from_true_anomaly,
                   cy_eccentric_anomaly_from_mean_anomaly_Newton1)
from .utils import ArrayProcessor

__all__ = ['mean_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_mean_anomaly',
           'true_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_true_anomaly']


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
    # TODO: in principle, this could be a decorator
    p = ArrayProcessor(E, e)
    E, e = p.prepare_arrays()
    return p.prepare_result(cy_mean_anomaly_from_eccentric_anomaly(E, e))


def eccentric_anomaly_from_mean_anomaly(M, e, tol=1E-10, maxiter=128):
    """
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
    method : str, optional
        The method to use for iterative root-finding for the eccentric anomaly.
        Options are: ``'Newton1'``.

    Returns
    -------
    E : numeric [radian]
        Eccentric anomaly.

    Issues
    ------
    - Magic numbers ``tol`` and ``maxiter``
    """

    p = ArrayProcessor(M, e)
    M, e = p.prepare_arrays()
    return p.prepare_result(cy_eccentric_anomaly_from_mean_anomaly_Newton1(
        M, e, tol, maxiter))


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
    f : numeric [radian]
        True anomaly.
    """
    p = ArrayProcessor(E, e)
    E, e = p.prepare_arrays()
    return p.prepare_result(cy_true_anomaly_from_eccentric_anomaly(E, e))


def eccentric_anomaly_from_true_anomaly(f, e):
    """
    Parameters
    ----------
    f : numeric, array_like [radian]
        True anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    E : numeric [radian]
        Eccentric anomaly.
    """
    p = ArrayProcessor(f, e)
    E, e = p.prepare_arrays()
    return p.prepare_result(cy_eccentric_anomaly_from_true_anomaly(E, e))

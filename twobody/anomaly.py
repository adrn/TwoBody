# Standard library
import inspect

# Third-party
import astropy.units as u
from astropy.utils.decorators import wraps

# Project
from .wrap import (cy_mean_anomaly_from_eccentric_anomaly,
                   cy_true_anomaly_from_eccentric_anomaly,
                   cy_eccentric_anomaly_from_true_anomaly,
                   cy_eccentric_anomaly_from_mean_anomaly_Newton1,
                   cy_eccentric_anomaly_from_mean_anomaly_Householder3)
from .utils import ArrayProcessor

__all__ = ['mean_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_mean_anomaly',
           'true_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_true_anomaly']


def anomaly_wrapper(func):
    """This is a decorator that processes the input (strips units, enforces that
    it is a 1D array), executes the cython function on the cleaned input, then
    returns the output angle in the same units and with the same shape as the
    input.
    """
    sig = inspect.signature(func)

    @wraps(func)
    def function_wrapper(ang, e, *args, **kwargs):
        in_unit = ang.unit
        p = ArrayProcessor(ang.to(u.radian).value, e)
        ang, e = p.prepare_arrays()
        res = p.prepare_result(func(ang, e, *args, **kwargs))
        return (res * u.radian).to(in_unit)

    return function_wrapper


@u.quantity_input(E=u.radian)
@anomaly_wrapper
def mean_anomaly_from_eccentric_anomaly(E, e):
    """
    Parameters
    ----------
    E : quantity_like [angle]
        Eccentric anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    M : numeric, array_like [radian]
        Mean anomaly.
    """
    return cy_mean_anomaly_from_eccentric_anomaly(E, e)


@u.quantity_input(M=u.radian)
@anomaly_wrapper
def eccentric_anomaly_from_mean_anomaly(M, e, tol=1E-10, maxiter=128,
                                        method='Newton1'):
    """
    Parameters
    ----------
    M : quantity_like [angle]
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
        Options are: ``'Newton1'`` and ``'Householder3'``.

    Returns
    -------
    E : numeric [radian]
        Eccentric anomaly.

    Issues
    ------
    - Magic numbers ``tol`` and ``maxiter``
    """
    func = eval('cy_eccentric_anomaly_from_mean_anomaly_{0}'.format(method))
    return func(M, e, tol, maxiter)


@u.quantity_input(E=u.radian)
@anomaly_wrapper
def true_anomaly_from_eccentric_anomaly(E, e):
    """
    Parameters
    ----------
    E : quantity_like [angle]
        Eccentric anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    f : numeric [radian]
        True anomaly.
    """
    return cy_true_anomaly_from_eccentric_anomaly(E, e)


@u.quantity_input(f=u.radian)
@anomaly_wrapper
def eccentric_anomaly_from_true_anomaly(f, e):
    """
    Parameters
    ----------
    f : quantity_like [angle]
        True anomaly.
    e : numeric, array_like
        Eccentricity.

    Returns
    -------
    E : numeric [radian]
        Eccentric anomaly.
    """
    return cy_eccentric_anomaly_from_true_anomaly(f, e)

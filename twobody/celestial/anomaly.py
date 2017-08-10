# Third-party
from astropy.utils.misc import check_broadcast
import numpy as np

from .wrap import (cy_mean_anomaly_from_eccentric_anomaly,
                   cy_eccentric_anomaly_from_mean_anomaly_Newton1,
                   cy_eccentric_anomaly_from_mean_anomaly_Householder3)

__all__ = ['mean_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_mean_anomaly',
           'true_anomaly_from_eccentric_anomaly',
           'eccentric_anomaly_from_true_anomaly',
           'd_eccentric_anomaly_d_mean_anomaly',
           'd_true_anomaly_d_eccentric_anomaly']

class _ArrayProcessor(object):

    def __init__(self, *arrs):
        self.arrs = [np.array(arr) for arr in arrs]

    def prepare_arrays(self):
        """Make sure input arrays are all C-contiguous and have same shape."""
        self.max_shape = None
        for arr in self.arrs:
            if self.max_shape is None:
                self.max_shape = arr.shape
            elif arr.shape > self.max_shape:
                self.max_shape = arr.shape

        orig_shapes = []
        arrs_1d = []
        for arr in self.arrs:
            orig_shapes.append(arr.shape)
            arr = np.broadcast_to(arr, self.max_shape).ravel()
            arrs_1d.append(np.ascontiguousarray(arr))

        if not check_broadcast(orig_shapes):
            raise ValueError("Shapes are not broadcastable: {0}"
                             .format(orig_shapes))

        return arrs_1d

    def prepare_result(self, res):
        return res.reshape(self.max_shape)

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
    p = _ArrayProcessor(E, e)
    E, e = p.prepare_arrays()
    return p.prepare_result(cy_mean_anomaly_from_eccentric_anomaly(E, e))

def eccentric_anomaly_from_mean_anomaly(M, e, tol=1E-10, maxiter=128,
                                        method='Newton1'):
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

    func_name = "cy_eccentric_anomaly_from_mean_anomaly_{0}".format(method)
    func = eval(func_name)

    p = _ArrayProcessor(M, e)
    M, e = p.prepare_arrays()
    return p.prepare_result(func(M, e, tol, maxiter))

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
    f = 2 * np.arctan2(np.sqrt(1+e) * np.sin(E/2),
                       np.sqrt(1-e) * np.cos(E/2))
    return f % (2*np.pi)

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
    E = np.arctan2(np.sqrt(1 - e**2) * np.sin(f),
                   e + np.cos(f))
    return E % (2*np.pi)

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

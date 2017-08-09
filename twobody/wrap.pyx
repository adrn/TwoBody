# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False

from __future__ import division, print_function

# Third-party
import numpy as np
cimport numpy as np
np.import_array()
import cython
cimport cython

# from libc.stdio cimport printf

cdef extern from "src/twobody.h":
    double c_mean_anomaly_from_eccentric_anomaly(double E, double e) nogil
    double c_eccentric_anomaly_from_mean_anomaly_Newton1(double M, double e,
                                                         double tol,
                                                         int maxiter) nogil

__all__ = ['_mean_anomaly_from_eccentric_anomaly',
           '_eccentric_anomaly_from_mean_anomaly_Newton1']

cpdef _mean_anomaly_from_eccentric_anomaly(
        np.ndarray[double, mode="c", ndim=1] E,
        np.ndarray[double, mode="c", ndim=1] e):
    """

    """
    cdef:
        double _E
        double _e
        int n
        int N = len(E)
        np.ndarray[double, mode="c", ndim=1] M = np.zeros(N)

    for n in range(N):
        M[n] = c_mean_anomaly_from_eccentric_anomaly(_E, _e)

    return M

cpdef _eccentric_anomaly_from_mean_anomaly_Newton1(
        np.ndarray[double, mode="c", ndim=1] M,
        np.ndarray[double, mode="c", ndim=1] e,
        double tol, int maxiter):
    """

    """
    cdef:
        double _M
        double _e
        int n
        int N = len(M)
        np.ndarray[double, mode="c", ndim=1] E = np.zeros(N)

    for n in range(N):
        E[n] = c_eccentric_anomaly_from_mean_anomaly_Newton1(_M, _e)

    return M


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
    double c_mean_anomaly_from_eccentric_anomaly(double E, double e)
    double c_eccentric_anomaly_from_mean_anomaly_Newton1(double M, double e,
                                                         double tol,
                                                         int maxiter)
    double c_eccentric_anomaly_from_mean_anomaly_Householder3(double M,
                                                              double e,
                                                              double tol,
                                                              int maxiter)

__all__ = ['cy_mean_anomaly_from_eccentric_anomaly',
           'cy_eccentric_anomaly_from_mean_anomaly_Newton1',
           'cy_eccentric_anomaly_from_mean_anomaly_Householder3']

cpdef cy_mean_anomaly_from_eccentric_anomaly(
        np.ndarray[double, mode="c", ndim=1] E,
        np.ndarray[double, mode="c", ndim=1] e):
    """

    """
    cdef:
        int n
        int N = len(E)
        np.ndarray[double, mode="c", ndim=1] M = np.zeros(N, dtype=np.float64)

    for n in range(N):
        M[n] = c_mean_anomaly_from_eccentric_anomaly(E[n], e[n])

    return M

cpdef cy_eccentric_anomaly_from_mean_anomaly_Newton1(
        np.ndarray[double, mode="c", ndim=1] M,
        np.ndarray[double, mode="c", ndim=1] e,
        double tol, int maxiter):
    """

    """
    cdef:
        int n
        int N = len(M)
        np.ndarray[double, mode="c", ndim=1] E = np.zeros(N, dtype=np.float64)

    for n in range(N):
        E[n] = c_eccentric_anomaly_from_mean_anomaly_Newton1(M[n], e[n],
                                                             tol, maxiter)

    return E

cpdef cy_eccentric_anomaly_from_mean_anomaly_Householder3(
        np.ndarray[double, mode="c", ndim=1] M,
        np.ndarray[double, mode="c", ndim=1] e,
        double tol, int maxiter):
    """

    """
    cdef:
        int n
        int N = len(M)
        np.ndarray[double, mode="c", ndim=1] E = np.zeros(N, dtype=np.float64)

    for n in range(N):
        E[n] = c_eccentric_anomaly_from_mean_anomaly_Householder3(M[n], e[n],
                                                                  tol, maxiter)

    return E


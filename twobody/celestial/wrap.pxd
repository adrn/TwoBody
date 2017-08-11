import numpy as np
cimport numpy as np
np.import_array()

cpdef cy_mean_anomaly_from_eccentric_anomaly(
    np.ndarray[double, mode="c", ndim=1] E,
    np.ndarray[double, mode="c", ndim=1] e)

cpdef cy_eccentric_anomaly_from_mean_anomaly_Newton1(
    np.ndarray[double, mode="c", ndim=1] M,
    np.ndarray[double, mode="c", ndim=1] e,
    double tol, int maxiter)

cpdef cy_eccentric_anomaly_from_mean_anomaly_Householder3(
    np.ndarray[double, mode="c", ndim=1] M,
    np.ndarray[double, mode="c", ndim=1] e,
    double tol, int maxiter)

cpdef cy_true_anomaly_from_eccentric_anomaly(
    np.ndarray[double, mode="c", ndim=1] E,
    np.ndarray[double, mode="c", ndim=1] e)

cpdef cy_eccentric_anomaly_from_true_anomaly(
    np.ndarray[double, mode="c", ndim=1] E,
    np.ndarray[double, mode="c", ndim=1] e)

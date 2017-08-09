#include <math.h>
#include "twobody.h"

double c_mean_anomaly_from_eccentric_anomaly(double E, double e) {
    /*
    Parameters
    ----------
    E : double [radian]
        Eccentric anomaly.
    e : double
        Eccentricity.

    Returns
    -------
    M : double [radian]
        Mean anomaly.
    */
    return E - e * sin(E);
}

double c_eccentric_anomaly_from_mean_anomaly_Newton1(double M, double e,
                                                     double tol, int maxiter) {
    /*
    tol = 1E-13
    maxiter = 128

    Parameters
    ----------
    Ms : double [radian]
        Mean anomaly.
    e : double
        Eccentricity.
    tol : double
        Numerical tolerance used in iteratively solving for eccentric anomaly.
    maxiter : int
        Maximum number of iterations when iteratively solving for eccentric
        anomaly.

    Returns
    -------
    E : double [radian]
        Eccentric anomaly.

    */

    double dM, E;

    // initialization
    E = M + e * sin(M);

    for (int i=0; i<maxiter; i++) {
        dM = M - c_mean_anomaly_from_eccentric_anomaly(E, e);
        E = E + dM / (1. - e * cos(E));

        if (abs(dM) < tol)
            break;
    }
}

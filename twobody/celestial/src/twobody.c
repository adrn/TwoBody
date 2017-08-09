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
    Parameters
    ----------
    M : double [radian]
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

        if (fabs(dM) < tol) {
            // printf("%d iterations\n", i);
            break;
        }

    }
    return E;
}

double c_eccentric_anomaly_from_mean_anomaly_Householder3(double M, double e,
                                                          double tol,
                                                          int maxiter) {
    /*
    Parameters
    ----------
    M : double [radian]
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

    Notes
    -----
    - uses 3rd-order Householder's method
    - follows exactly https://en.wikipedia.org/wiki/Householder%27s_method

    */

    double dM, E;
    double f, fp, fpp, fppp, h;
    double ecosE, esinE;

    // initialization
    E = M + e * sin(M);

    for (int i=0; i<maxiter; i++) {
        ecosE = e * cos(E);
        esinE = e * sin(E);

        f = E - esinE - M; // function
        fp = 1 - ecosE;    // first derivative
        fpp = esinE;       // second
        fppp = ecosE;      // third
        h = -f / fp;       // Newton step

        E = E + h * (1. + (fpp / (2. * fp)) * h)
            / (1. + (fpp / fp) * h + (fppp / (6. * fp)) * h * h);

        if (fabs(f) < tol) {
            // printf("%d iterations\n", i);
            break;
        }

    }
    return E;
}

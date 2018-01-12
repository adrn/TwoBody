#include <math.h>
#include "twobody.h"

double mod_angle(double x){
    x = fmod(x, 2*M_PI);
    if (x < 0)
        x += 2*M_PI;
    return x;
}

// Mean anomaly <--> Eccentric anomaly

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

// True anomaly <--> Eccentric anomaly

double c_true_anomaly_from_eccentric_anomaly(double E, double e) {
    /*
    Parameters
    ----------
    E : double [radian]
        Eccentric anomaly.
    e : double
        Eccentricity.

    Returns
    -------
    f : double [radian]
        True anomaly.
    */
    double f;

    f = 2 * atan2(sqrt(1+e) * sin(E/2),
                  sqrt(1-e) * cos(E/2));

    return mod_angle(f);
}

double c_eccentric_anomaly_from_true_anomaly(double f, double e) {
    /*
    Parameters
    ----------
    f : double [radian]
        True anomaly.
    e : double
        Eccentricity.

    Returns
    -------
    E : double [radian]
        Eccentric anomaly.

    */
    double E;
    E = atan2(sqrt(1 - e*e) * sin(f), e + cos(f));
    return mod_angle(E);
}

void c_rv_from_elements(double *t, double *rv, int N_t,
                        double P, double K, double e, double omega,
                        double M0, double t0, double tol, int maxiter) {
    /* Compute the RV of the primary w.r.t. the system barycenter */
    double M, E, f;
    for (int n=0; n < N_t; n++) {
        M = 2 * M_PI * (t[n] - t0) / P - M0;
        E = c_eccentric_anomaly_from_mean_anomaly_Newton1(M, e,
                                                          tol, maxiter);
        f = c_true_anomaly_from_eccentric_anomaly(E, e);
        rv[n] = K * (cos(omega + f) + e * cos(omega));
    }
}

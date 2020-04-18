#ifndef M_PI
    #define M_PI 3.14159265358979323846264338327950288
#endif

extern double c_mean_anomaly_from_eccentric_anomaly(double E, double e);
extern double c_eccentric_anomaly_from_mean_anomaly_Newton1(double M, double e,
                                                            double tol,
                                                            int maxiter);
extern double c_eccentric_anomaly_from_mean_anomaly_Householder3(double M,
                                                                 double e,
                                                                 double tol,
                                                                 int maxiter);

extern double c_true_anomaly_from_eccentric_anomaly(double E, double e);
extern double c_eccentric_anomaly_from_true_anomaly(double f, double e);

extern void c_rv_from_elements(double *t, double *rv, int N_t,
                               double P, double K, double e, double omega,
                               double M0, double t0, double tol, int maxiter);

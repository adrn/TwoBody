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

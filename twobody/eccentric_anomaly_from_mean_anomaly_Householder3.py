def eccentric_anomaly_from_mean_anomaly_Householder3(Ms, e, tol=1E-13, maxiter=128):
    """
    Parameters
    ----------
    Ms : numeric, array_like [radian]
        Mean anomaly.
    e : numeric
        Eccentricity.
    tol : numeric, optional
        Numerical tolerance used in iteratively solving for eccentric anomaly.
    maxiter : int, optional
        Maximum number of iterations when iteratively solving for
        eccentric anomaly.

    Returns
    -------
    Es : numeric [radian]
        Eccentric anomaly.

    Notes
    -----
    - uses 3rd-order Householder's method
    - follows exactly https://en.wikipedia.org/wiki/Householder%27s_method

    Issues
    ------
    - magic numbers ``tol`` and ``maxiter``
    """
    Ms = np.atleast_1d(Ms)

    if Ms.ndim > 1:
        raise ValueError("Input must have <= 1 dim.")

    Es = Ms + e * np.sin(Ms)   # first guess
    bad = Es < np.Inf          # list of bools
    for _ in range(maxiter):
        ecosEs, esinEs = e * np.cos(Es[bad]), e * np.sin(Es[bad])  # pre-compute
        f = Es[bad] - esinEs - Ms[bad]  # function
        fp = 1. - ecosEs       # first derivative
        fpp = esinEs           # second derivative
        fppp = ecosEs          # third
        h = - f / fp           # Newton step
        Es[bad] = Es[bad] + h * (1. + (fpp / (2. * fp)) * h) \
            / (1. + (fpp / fp) * h + (fppp / (6. * fp)) * h * h)

        bad[bad] = f > tol  # check my indexing shih
        if not np.any(bad):
            break
    else:
        warnings.warn("eccentric_anomaly_from_mean_anomaly() reached maximum "
                      "number of iterations ({})".format(maxiter), RuntimeWarning)
    return Es


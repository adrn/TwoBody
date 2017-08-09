# Third-party
import numpy as np
import pytest

# Project
from .. import anomaly
from ..anomaly import (mean_anomaly_from_eccentric_anomaly,
                       eccentric_anomaly_from_mean_anomaly)

N = 16

def test_mean_anomaly_from_eccentric_anomaly():
    E_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1, N)
    Es, es = [x.ravel() for x in np.meshgrid(E_grid, e_grid)]

    # call with float
    for E, e in zip(Es, es):
        M = mean_anomaly_from_eccentric_anomaly(E, e)

    # call with arrays
    M = mean_anomaly_from_eccentric_anomaly(Es, es)

@pytest.mark.parametrize("method", ["Newton1"])
def test_eccentric_anomaly_from_mean_anomaly(method):
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1, N)
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    # call with float
    for M, e in zip(Ms, es):
        E = eccentric_anomaly_from_mean_anomaly(M, e, method=method)

    # call with arrays
    E = eccentric_anomaly_from_mean_anomaly(Ms, es)

@pytest.mark.parametrize("method", ["Newton1"])
def test_anomaly_roundtrip(method):
    func = getattr(anomaly,
                   "eccentric_anomaly_from_mean_anomaly_{0}".format(method))
    kw = dict(tol=1E-10, maxiter=128)

    # M -> E -> M
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-9, N) # MAGIC NUMBER: 1E-9
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    for M, e in zip(Ms, es):
        E = func(M, e, **kw)
        M2 = mean_anomaly_from_eccentric_anomaly(E, e)
        assert np.allclose(M, M2, atol=1E-14)

    E = func(Ms, es, **kw)
    M2 = mean_anomaly_from_eccentric_anomaly(E, es)
    assert np.allclose(Ms, M2, atol=1E-14)

    # E -> M -> E
    E_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-9, N) # MAGIC NUMBER: 1E-9
    Es, es = [x.ravel() for x in np.meshgrid(E_grid, e_grid)]

    for E, e in zip(Es, es):
        M = mean_anomaly_from_eccentric_anomaly(E, e)
        E2 = func(M, e, **kw)
        assert np.allclose(E, E2, atol=1E-14)

    M = mean_anomaly_from_eccentric_anomaly(Es, es)
    E2 = func(M, es, **kw)
    assert np.allclose(Es, E2, atol=1E-14)


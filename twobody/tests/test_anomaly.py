# Third-party
from astropy.tests.helper import quantity_allclose
import astropy.units as u
import numpy as np
import pytest

# Project
from ..anomaly import (mean_anomaly_from_eccentric_anomaly,
                       eccentric_anomaly_from_mean_anomaly,
                       true_anomaly_from_eccentric_anomaly,
                       eccentric_anomaly_from_true_anomaly)

N = 16


def test_mean_anomaly_from_eccentric_anomaly():
    E_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1, N)
    Es, es = [x.ravel() for x in np.meshgrid(E_grid, e_grid)]

    # call with float
    for E, e in zip(Es, es):
        M = mean_anomaly_from_eccentric_anomaly(E*u.rad, e)
        assert M.shape == ()

    # call with arrays
    M = mean_anomaly_from_eccentric_anomaly(Es*u.rad, es)
    assert M.shape == Es.shape


@pytest.mark.parametrize("method", ["Newton1", "Householder3"])
def test_eccentric_anomaly_from_mean_anomaly(method):
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1, N)
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    # call with float
    for M, e in zip(Ms, es):
        E = eccentric_anomaly_from_mean_anomaly(M*u.rad, e, method=method)
        assert E.shape == ()

    # call with arrays
    E = eccentric_anomaly_from_mean_anomaly(Ms*u.rad, es)
    assert E.shape == Ms.shape


@pytest.mark.parametrize("method", ["Newton1", "Householder3"])
def test_anomaly_roundtrip(method):
    kw = dict(tol=1E-10, maxiter=128, method=method)

    # M -> E -> M
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-3, N)  # MAGIC NUMBER: 1E-3
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    for M, e in zip(Ms, es):
        E = eccentric_anomaly_from_mean_anomaly(M*u.rad, e, **kw)
        M2 = mean_anomaly_from_eccentric_anomaly(E, e)
        assert quantity_allclose(M*u.rad, M2, atol=1E-14*u.rad)

    E = eccentric_anomaly_from_mean_anomaly(Ms*u.rad, es, **kw)
    M2 = mean_anomaly_from_eccentric_anomaly(E, es)
    assert quantity_allclose(Ms*u.rad, M2, atol=1E-14*u.rad)

    # E -> M -> E
    E_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-3, N)  # MAGIC NUMBER: 1E-3
    Es, es = [x.ravel() for x in np.meshgrid(E_grid, e_grid)]

    for E, e in zip(Es, es):
        M = mean_anomaly_from_eccentric_anomaly(E*u.rad, e)
        E2 = eccentric_anomaly_from_mean_anomaly(M, e, **kw)
        assert quantity_allclose(E*u.rad, E2, atol=1E-14*u.rad)

    M = mean_anomaly_from_eccentric_anomaly(Es*u.rad, es)
    E2 = eccentric_anomaly_from_mean_anomaly(M, es, **kw)
    assert quantity_allclose(Es*u.rad, E2, atol=1E-14*u.rad)


def test_true_anomaly_roundtrip():
    # f -> E -> f
    f_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-3, N)  # MAGIC NUMBER: 1E-3
    fs, es = [x.ravel() for x in np.meshgrid(f_grid, e_grid)]

    E = eccentric_anomaly_from_true_anomaly(fs*u.rad, es)
    f2 = true_anomaly_from_eccentric_anomaly(E, es)
    assert np.allclose(np.cos(fs), np.cos(f2), atol=1E-9)

    # E -> f -> E
    E_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-3, N)  # MAGIC NUMBER: 1E-3
    Es, es = [x.ravel() for x in np.meshgrid(E_grid, e_grid)]

    f = true_anomaly_from_eccentric_anomaly(Es*u.rad, es)
    E2 = eccentric_anomaly_from_true_anomaly(f, es)
    assert np.allclose(np.cos(Es), np.cos(E2), atol=1E-9)

# Standard library
import os
import sys

# Third-party
import numpy as np

# Project
from ..core import (mean_anomaly_from_eccentric_anomaly,
                    eccentric_anomaly_from_mean_anomaly_Newton1)

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

def test_eccentric_anomaly_from_mean_anomaly_Newton1():
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1, N)
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    # call with float
    for M, e in zip(Ms, es):
        E = mean_anomaly_from_eccentric_anomaly(M, e)

    # call with arrays
    E = mean_anomaly_from_eccentric_anomaly(Ms, es)

def test_anomaly_roundtrip():
    # M -> E -> M
    M_grid = np.linspace(0, 2*np.pi, N)
    e_grid = np.linspace(0, 1 - 1E-9, N) # MAGIC NUMBER: 1E-9
    Ms, es = [x.ravel() for x in np.meshgrid(M_grid, e_grid)]

    for M, e in zip(Ms, es):
        E = eccentric_anomaly_from_mean_anomaly_Newton1(M, e, tol=1E-10)
        M2 = mean_anomaly_from_eccentric_anomaly(E, e)
        assert np.allclose(M, M2, atol=1E-14)

    E = eccentric_anomaly_from_mean_anomaly_Newton1(Ms, es, tol=1E-10)
    M2 = mean_anomaly_from_eccentric_anomaly(E, es)
    assert np.allclose(Ms, M2, atol=1E-14)


# Third-party
import astropy.units as u
from astropy.time import Time
import numpy as np
import pytest

from ..trends import (VelocityTrend1, VelocityTrend2, VelocityTrend3)

coeffs = [10.*u.km/u.s, 1.*u.km/u.s/u.day, 1.*u.km/u.s/u.day**2]
coeffs_no_u = [1., 1, 1]
bad_coeffs = [1.*u.km/u.s/u.day, 1.*u.km/u.s/u.day**2, 1.*u.km/u.s/u.day**3]

@pytest.mark.parametrize("cls", [VelocityTrend1,
                                 VelocityTrend2,
                                 VelocityTrend3])
def test_subclasses(cls):
    # initialization
    n_pars = len(cls.parameters)
    trend = cls(*coeffs[:n_pars])

    # invalid units
    with pytest.raises(ValueError):
        cls(*coeffs_no_u[:n_pars])

    with pytest.raises(u.UnitsError):
        cls(*bad_coeffs[:n_pars])

    trend = cls(*coeffs[:n_pars])
    t = np.random.uniform(15., 23., 128)
    res1 = trend(t)
    res2 = trend(t*u.day)
    assert np.allclose(res1.value, res2.value)
    assert res1.unit == u.km/u.s
    assert res2.unit == u.km/u.s

    with pytest.raises(u.UnitsError):
        trend(np.random.uniform(15., 23., 128)*u.km)

    # Now with t0
    trend = cls(*coeffs[:n_pars], t0=57423.)
    res1 = trend(t*u.day)

    trend = cls(*coeffs[:n_pars], t0=Time(57423, format='mjd', scale='utc'))
    res2 = trend(t*u.day)
    assert np.allclose(res1.value, res2.value)

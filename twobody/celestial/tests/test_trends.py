# Third-party
import astropy.units as u
from astropy.time import Time
import numpy as np
import pytest

from ..trends import PolynomialVelocityTrend

def test_polynomial():
    # initialization
    trend = PolynomialVelocityTrend(10.*u.km/u.s,
                                    1.*u.km/u.s/u.day,
                                    1.*u.km/u.s/u.day**2)

    trend = PolynomialVelocityTrend()
    with pytest.raises(ValueError):
        trend(0.)

    # invalid units
    with pytest.raises(ValueError):
        PolynomialVelocityTrend(10., 1.*u.km/u.s/u.day)

    with pytest.raises(u.UnitsError):
        PolynomialVelocityTrend(10.*u.km/u.s, 1.*u.km/u.s)

    trend = PolynomialVelocityTrend(10.*u.km/u.s, 1.*u.km/u.s/u.day)
    t = np.random.uniform(15., 23., 128)
    res1 = trend(t)
    res2 = trend(t*u.day)
    assert np.allclose(res1.value, res2.value)
    assert res1.unit == u.km/u.s
    assert res2.unit == u.km/u.s

    with pytest.raises(u.UnitsError):
        trend(np.random.uniform(15., 23., 128)*u.km)

    # Now with t0
    trend = PolynomialVelocityTrend(10.*u.km/u.s,
                                    1.*u.km/u.s/u.day,
                                    1.*u.km/u.s/u.day**2,
                                    t0=57423.)
    res1 = trend(t*u.day)

    trend = PolynomialVelocityTrend(10.*u.km/u.s,
                                    1.*u.km/u.s/u.day,
                                    1.*u.km/u.s/u.day**2,
                                    t0=Time(57423, format='mjd', scale='utc'))
    res2 = trend(t*u.day)
    assert np.allclose(res1.value, res2.value)

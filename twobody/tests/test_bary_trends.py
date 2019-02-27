# Third-party
import astropy.units as u
from astropy.time import Time
import numpy as np
import pytest

from ..bary_trends import PolynomialRVTrend

def test_polynomial():
    # initialization
    trend = PolynomialRVTrend([10.*u.km/u.s,
                               1.*u.km/u.s/u.day,
                               1.*u.km/u.s/u.day**2])

    trend = PolynomialRVTrend()
    assert trend(0.) == 0.
    assert trend([0.]) == np.array([0.])
    assert np.all(trend(np.arange(10)) == 0.)

    # invalid units
    with pytest.raises(ValueError):
        PolynomialRVTrend([10., 1.*u.km/u.s/u.day])

    with pytest.raises(u.UnitsError):
        PolynomialRVTrend([10.*u.km/u.s, 1.*u.km/u.s])

    trend = PolynomialRVTrend([10.*u.km/u.s, 1.*u.km/u.s/u.day])
    t = Time(55555. + np.sort(np.random.uniform(15., 23., 128)),
             format='mjd')
    res1 = trend(t)
    res2 = trend(t.tcb.mjd)
    assert np.allclose(res1.value, res2.value)
    assert res1.unit == u.km/u.s
    assert res2.unit == u.km/u.s

    trend = PolynomialRVTrend([10.*u.km/u.s, 1.*u.km/u.s/u.day],
                              t0=Time(55555., format='mjd'))
    res = trend(t.tcb.mjd)

    t = (Time(55555., format='mjd') +
         np.sort(np.random.uniform(15., 23., 128))*u.day)
    res = trend(t)

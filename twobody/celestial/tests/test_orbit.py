# Third-party
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import quantity_allclose
import numpy as np
import pytest

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except:
    HAS_MPL = False

# Package
from ..orbit import RVOrbit
from ..trends import PolynomialVelocityTrend

# default orbit parameters
default_pars = dict(P=30.*u.day, K=100.*u.m/u.s, ecc=0.11239,
                    omega=0.*u.radian, phi0=0.25524*u.radian)

def test_init_rvorbit():

    # without trend
    orbit = RVOrbit(**default_pars)

    t = np.random.uniform(55612., 55792, 1024)
    t.sort()
    rv = orbit.generate_rv_curve(t)
    rv2 = orbit(t)
    assert quantity_allclose(rv, rv2)

    # with trend - check t0
    pars = dict(P=1.*u.day, K=100.*u.m/u.s, ecc=0.,
                omega=0.*u.radian, phi0=90*u.degree)
    trend = PolynomialVelocityTrend(100.*u.km/u.s, 1*u.km/u.s/u.day, t0=55831.)
    orbit1 = RVOrbit(trend=trend, **pars)
    assert quantity_allclose(orbit1(55831.), 100 * u.km/u.s)

    # get pericenter time - check something?!
    orbit = RVOrbit(**pars)
    mjd_grid = [57975.15, 57975.95, 57978.55]
    expected = [57975.25, 57976.25, 57978.25]
    mjd0 = [orbit.t0(x).tcb.mjd for x in mjd_grid]
    assert np.allclose(mjd0, expected)

@pytest.mark.skipif(not HAS_MPL, reason="matplotlib not installed")
def test_plotting():
    orbit = RVOrbit(**default_pars)

    mjd = np.linspace(56823.123, 57293.2345, 1024) # random MJD's
    t = Time(mjd, format='mjd', scale='utc')

    # plotting
    for _t in [mjd, t]:
        orbit.plot(_t)
        orbit.plot(_t, t_kwargs=dict(format='jd', scale='utc'))

    fig,ax = plt.subplots(1,1)
    orbit.plot(t=t, ax=ax, plot_kwargs=dict(color='#de2d26'))

    fig,ax = plt.subplots(1,1)
    orbit.plot(t=mjd, ax=ax)

    plt.close('all')

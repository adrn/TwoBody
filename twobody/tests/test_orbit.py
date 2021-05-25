# Third-party
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.tests.helper import quantity_allclose
import numpy as np
import pytest

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# Package
from ..orbit import KeplerOrbit
from ..barycenter import Barycenter
from ..bary_trends import PolynomialRVTrend

def get_random_orbit_pars(rnd=None, barycenter=None):
    if rnd is None:
        rnd = np.random.RandomState()

    P = 10**rnd.uniform(0, 3)*u.day
    K = np.abs(rnd.normal(0, 100))*u.km/u.s
    e = rnd.uniform(0, 1)
    i = rnd.uniform(0, 180)*u.deg
    n = 2*np.pi/P
    a = K * np.sqrt(1-e**2) / (n * np.sin(i))
    return dict(P=P, a=a, e=e, i=i,
                omega=rnd.uniform(0, 360)*u.deg,
                Omega=rnd.uniform(0, 360)*u.deg,
                M0=rnd.uniform(0, 360)*u.deg,
                barycenter=barycenter)


def test_radial_velocity():
    # Compare the short-circuit RV to the full calculation

    rnd = np.random.RandomState(seed=42)
    n_iter = 32

    bc_kw = dict(ra=210.413231*u.deg,
                 dec=-61.3412*u.deg,
                 distance=101*u.pc)

    # at rest w.r.t. sun
    bc1 = Barycenter(origin=coord.ICRS(**bc_kw,
                                       pm_ra_cosdec=0*u.mas/u.yr,
                                       pm_dec=0.*u.mas/u.yr,
                                       radial_velocity=0*u.km/u.s),
                     t0=Time('J2018.0'))

    bc2 = Barycenter(origin=coord.ICRS(**bc_kw,
                                       pm_ra_cosdec=1000*u.mas/u.yr,
                                       pm_dec=1000.*u.mas/u.yr,
                                       radial_velocity=-260*u.km/u.s),
                     t0=Time('J2018.0'))

    times = Time('J2018.0') + np.linspace(-100, 100, 128) * u.day

    # Even when barycenter is moving, rv1 and rv2 should be the same
    for n in range(n_iter):
        orb = KeplerOrbit(**get_random_orbit_pars(rnd, bc2))

        rv1 = orb.radial_velocity(times)
        _rp = orb.reference_plane(times)
        rv2 = -_rp.cartesian.differentials['s'].d_z + bc2.origin.radial_velocity

        assert quantity_allclose(rv1, rv2, rtol=0, atol=1E-3*u.m/u.s)

    # All RV's should be equivalent
    for n in range(n_iter):
        orb = KeplerOrbit(**get_random_orbit_pars(rnd, bc1))

        rv1 = orb.radial_velocity(times)
        _rp = orb.reference_plane(times)
        rv2 = -_rp.cartesian.differentials['s'].d_z
        rv3 = orb.icrs(times).radial_velocity

        assert quantity_allclose(rv1, rv2, rtol=0, atol=1E-4*u.m/u.s)
        # We expect true differences because of projection effects
        assert quantity_allclose(rv1, rv3, rtol=0, atol=1*u.km/u.s)


def test_velocity_trend():
    rnd = np.random.RandomState(seed=42)

    coeffs = [0*u.km/u.s,
              0*u.km/u.s/u.day,
              1e-2*u.km/u.s/u.day**2]
    trend = PolynomialRVTrend(coeffs,
                              t0=Time('J2018.0'))
    pars = get_random_orbit_pars(rnd, trend)
    pars['a'] = 0*u.au
    orbit = KeplerOrbit(**pars)

    t = Time('J2018.0') + np.linspace(-100, 100, 256)*u.day
    rv = orbit.radial_velocity(t)
    assert quantity_allclose(rv[0], rv[-1])
    assert quantity_allclose(rv[0], 100.*u.km/u.s)

    # OK THIS SUCKS (copy pasta)
    coeffs = [0*u.km/u.s,
              1*u.km/u.s/u.day,
              0*u.km/u.s/u.day**2]
    trend = PolynomialRVTrend(coeffs,
                              t0=Time('J2018.0'))
    pars = get_random_orbit_pars(rnd, trend)
    pars['a'] = 0*u.au
    orbit = KeplerOrbit(**pars)

    t = Time('J2018.0') + np.linspace(-100, 100, 256)*u.day
    rv = orbit.radial_velocity(t)
    assert quantity_allclose(rv[0], -rv[-1])
    assert quantity_allclose(rv[0], -100.*u.km/u.s)


@pytest.mark.skipif(not HAS_MPL, reason="matplotlib not installed")
def test_plotting():
    orbit = KeplerOrbit(P=100*u.day, a=1.*u.au, e=0.1,
                        omega=0*u.deg, i=41*u.deg, Omega=0*u.deg)

    mjd = np.linspace(56823.123, 57293.2345, 1024)  # random MJD's
    t = Time(mjd, format='mjd', scale='utc')

    # plotting
    for _t in [mjd, t]:
        orbit.plot_rv(_t)
        orbit.plot_rv(_t, t_kwargs=dict(format='jd', scale='utc'))

    fig,ax = plt.subplots(1,1)
    orbit.plot_rv(time=t, ax=ax, plot_kwargs=dict(color='#de2d26'))

    fig,ax = plt.subplots(1,1)
    orbit.plot_rv(time=mjd, ax=ax)

    plt.close('all')


def test_nan():

    mjd = np.linspace(56123.123, 57293.2345, 1024)  # random MJD's
    t = Time(mjd, format='mjd', scale='utc')

    orb = KeplerOrbit(P=1.5*u.year, e=0.67, M0=0*u.deg, omega=17.14*u.deg,
                      i=np.nan*u.deg, Omega=np.nan*u.deg, a=np.nan*u.au)
    rv = orb.radial_velocity(t)
    urv = orb.unscaled_radial_velocity(t)

    assert np.all(np.isnan(rv.value))
    assert np.all(np.isfinite(urv))

    orb = KeplerOrbit(P=1.5*u.year, e=0.67, M0=0*u.deg, omega=17.14*u.deg)
    with pytest.raises(ValueError):
        rv = orb.radial_velocity(t)
    urv = orb.unscaled_radial_velocity(t)

    assert np.all(np.isnan(rv.value))
    assert np.all(np.isfinite(urv))

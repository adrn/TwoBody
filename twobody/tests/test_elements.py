# Third-party
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import catch_warnings
import numpy as np
import pytest

# Project
from ..elements import KeplerElements, TwoBodyKeplerElements
from ..orbit import KeplerOrbit
from ..units import UnitSystem


def test_kepler():

    # Valid initialization:
    kw1 = dict(P=10*u.day, e=0.5,
               omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)
    kw2 = dict(a=0.1*u.au)
    kw3 = dict(M0=40*u.deg, t0=Time('J2015.0'))
    kw4 = dict(units=UnitSystem(u.m, u.s, u.rad, u.kg))
    kw5 = dict(K=10*u.km/u.s)

    for kw in [kw1, {**kw1, **kw2}, {**kw1, **kw2, **kw3}, {**kw1, **kw3},
               {**kw1, **kw4}, {**kw1, **kw2, **kw4},
               {**kw1, **kw2, **kw3, **kw4},
               {**kw1, **kw5}]:  # not meant to be exhaustive
        elems = KeplerElements(**kw)

        # Make sure this also works as arguments to KeplerOrbit
        KeplerOrbit(elems, elements_type='kepler')
        KeplerOrbit(**kw, elements_type='kepler')

    # Invalid initialization:
    for name in ['P', 'omega']:
        kw = kw1.copy()
        kw.pop(name)
        with pytest.raises(ValueError):
            KeplerElements(**kw)

    # check they get set to nan:
    kw = kw1.copy()
    kw.pop('Omega')
    elem = KeplerElements(**kw)
    assert np.isnan(elem.Omega)

    kw = kw1.copy()
    kw.pop('i')
    elem = KeplerElements(**kw)
    assert np.isnan(elem.i)

    # Make sure attributes exist and have correct units:
    elems = KeplerElements(P=10*u.day, e=0.5, a=0.1*u.au,
                           omega=10*u.deg, i=20*u.deg, Omega=30*u.deg,
                           units=UnitSystem(u.m, u.s, u.rad, u.kg, u.m/u.s))

    assert elems.P.unit == u.s
    assert elems.K.unit == u.m/u.s
    assert elems.m_f.unit == u.kg

    elems = KeplerElements(P=10*u.day, e=0.5, K=10*u.km/u.s,
                           omega=10*u.deg, i=20*u.deg, Omega=30*u.deg,
                           units=UnitSystem(u.m, u.s, u.rad, u.kg, u.m/u.s))

    assert elems.P.unit == u.s
    assert elems.K.unit == u.m/u.s
    assert elems.K == 10*u.km/u.s
    assert elems.m_f.unit == u.kg

    # Expected failures:
    bad_kws = [('P', -10*u.day), ('a', -10*u.au), ('e', 1.5), ('i', -10*u.deg)]
    for k, v in bad_kws:
        kw = kw1.copy()
        kw[k] = v
        with pytest.raises(ValueError):
            KeplerElements(**kw)

    with pytest.raises(ValueError):
        elems = KeplerElements(P=10*u.day, e=0.5, K=10*u.km/u.s, a=0.1*u.au,
                               omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)

    # Check warning
    with catch_warnings() as w:
        KeplerElements(P=1*u.day, a=5*u.au, e=0.,
                       omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)
    assert len(w) > 0


def test_twobodykepler():

    TwoBodyKeplerElements(P=10*u.day, e=0.5,
                          m1=1*u.Msun, m2=5*u.Msun,
                          omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)
    TwoBodyKeplerElements(a=10*u.au, e=0.5,
                          m1=1*u.Msun, m2=5*u.Msun,
                          omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)

    # check attributes
    elems = TwoBodyKeplerElements(P=10*u.day, e=0.5,
                                  m1=5*u.Msun, m2=1*u.Msun,
                                  omega=10*u.deg, i=20*u.deg, Omega=30*u.deg)

    with pytest.raises(AttributeError):
        elems.K

    with pytest.raises(AttributeError):
        elems.m_f

    elems.get_body('1')
    elems.get_body('2')

    elems.get_body(1)
    elems.get_body(2)

    with pytest.raises(ValueError):
        elems.get_body('derp')

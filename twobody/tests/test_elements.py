# Third-party
from astropy.time import Time
import astropy.units as u
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

    for kw in [kw1, {**kw1, **kw2}, {**kw1, **kw2, **kw3}, {**kw1, **kw3},
               {**kw1, **kw4}, {**kw1, **kw2, **kw4},
               {**kw1, **kw2, **kw3, **kw4}]: # not meant to be fully exhaustive
        elems = KeplerElements(**kw)

        # Make sure this also works as arguments to KeplerOrbit
        KeplerOrbit(elems, elements_type='kepler')
        KeplerOrbit(**kw, elements_type='kepler')

    # Invalid initialization:
    for name in ['P', 'i', 'omega', 'Omega']:
        kw = kw1.copy()
        kw.pop(name)
        with pytest.raises(ValueError):
            KeplerElements(**kw)

    # Make sure attributes exist and have correct units:
    elems = KeplerElements(P=10*u.day, e=0.5, a=0.1*u.au,
                           omega=10*u.deg, i=20*u.deg, Omega=30*u.deg,
                           units=UnitSystem(u.m, u.s, u.rad, u.kg, u.m/u.s))

    assert elems.P.unit == u.s
    assert elems.K.unit == u.m/u.s
    assert elems.m_f.unit == u.kg

    # Expected failures:
    bad_kws = [('P',-10*u.day), ('a',-10*u.au), ('e',1.5), ('i',-10*u.deg)]
    for k,v in bad_kws:
        kw = kw1.copy()
        kw[k] = v
        with pytest.raises(ValueError):
            KeplerElements(**kw)

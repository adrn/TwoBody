# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy.tests.helper import quantity_allclose
import numpy as np
import pytest

# Project
from ..reference_plane import ReferencePlaneFrame


def test_sanity():

    rep1 = coord.CartesianRepresentation(x=[0, 1.],
                                         y=0,
                                         z=0, unit=u.pc)

    rep2 = coord.CartesianRepresentation(x=0,
                                         y=[0, 1.],
                                         z=0, unit=u.pc)

    rep3 = coord.CartesianRepresentation(x=0,
                                         y=0,
                                         z=[0, 1.], unit=u.pc)

    # Try for many origins:
    rnd = np.random.RandomState(seed=42)

    for _ in range(128):
        origin = coord.ICRS(ra=rnd.uniform(0, 360)*u.deg,
                            dec=rnd.uniform(-90, 90)*u.deg,
                            distance=rnd.uniform(10, 100)*u.pc)

        ref_c = ReferencePlaneFrame(rep1, origin=origin)
        icrs = ref_c.transform_to(coord.ICRS)
        assert icrs.dec[1] > icrs.dec[0]
        assert quantity_allclose(icrs.ra[0], icrs.ra[1])

        ref_c = ReferencePlaneFrame(rep2, origin=origin)
        icrs = ref_c.transform_to(coord.ICRS)
        assert icrs.ra[1] > icrs.ra[0]

        ref_c = ReferencePlaneFrame(rep3, origin=origin)
        icrs = ref_c.transform_to(coord.ICRS)
        assert icrs.distance[0] > icrs.distance[1]
        assert quantity_allclose(icrs.ra[0], icrs.ra[1])
        assert quantity_allclose(icrs.dec[0], icrs.dec[1])


def test_transform():
    rep = coord.CartesianRepresentation(x=1, y=2, z=3, unit=u.pc)
    ref_c1 = ReferencePlaneFrame(rep)

    with pytest.raises(ValueError):
        ref_c1.transform_to(coord.ICRS)

    with pytest.raises(ValueError):
        ref_c2 = ReferencePlaneFrame(rep, origin=coord.Galactic())

# Third-party
import astropy.coordinates as coord
import astropy.units as u
import pytest
from astropy.time import Time

# Project
from ..barycenter import Barycenter


def test_init():
    b = Barycenter(
        origin=coord.SkyCoord(ra=150 * u.deg, dec=-11 * u.deg, distance=100 * u.pc)
    )

    b = Barycenter(
        origin=coord.SkyCoord(ra=150 * u.deg, dec=-11 * u.deg, distance=100 * u.pc),
        t0=59123.1213,
    )

    b = Barycenter(
        origin=coord.SkyCoord(ra=150 * u.deg, dec=-11 * u.deg, distance=100 * u.pc),
        t0=Time(59123.1213, format="mjd"),
    )

    c = coord.SkyCoord(ra=150 * u.deg, dec=-11 * u.deg, distance=100 * u.pc)
    with pytest.raises(ValueError):
        Barycenter(origin=c, radial_velocity=100 * u.km / u.s)

    with pytest.raises(ValueError):
        Barycenter()

    with pytest.raises(TypeError):
        Barycenter(origin=Time("J2000"))

    with pytest.raises(ValueError):
        Barycenter(
            origin=coord.SkyCoord(
                ra=[150, 100] * u.deg, dec=[-11, 10] * u.deg, distance=100 * u.pc
            )
        )

    with pytest.warns(UserWarning, match="Barycenter") as w:
        Barycenter(origin=coord.SkyCoord(ra=150 * u.deg, dec=-11 * u.deg))

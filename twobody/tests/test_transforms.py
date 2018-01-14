# Third-party
import astropy.units as u
from astropy.tests.helper import quantity_allclose
import numpy as np

# Project
from ..transforms import a_P_to_m, a_m_to_P, P_m_to_a


def test_roundtrip_Pma():

    rnd = np.random.RandomState(seed=42)

    for i in range(128):
        P = 10 ** rnd.uniform(0, 3) * u.day
        m = np.abs(rnd.normal(1., 0.5)) * u.Msun
        a = P_m_to_a(P, m)

        m2 = a_P_to_m(a, P)
        P2 = a_m_to_P(a, m2)
        a2 = P_m_to_a(P2, m2)

        assert quantity_allclose(a, a2)
        assert quantity_allclose(m, m2)
        assert quantity_allclose(P, P2)

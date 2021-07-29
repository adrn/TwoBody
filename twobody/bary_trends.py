import astropy.units as u
from astropy.time import Time
import numpy as np

__all__ = ['RVTrend', 'PolynomialRVTrend']


class RVTrend:
    pass


class PolynomialRVTrend(RVTrend):

    def __init__(self, coeffs=None, t0=None):
        """Specify a polynomial trend in radial velocities by specifying the
        polynomial coefficients and optionally a reference time.

        Note that the coefficients are expected in the reverse order of those
        expected by numpy's ``poly`` functionality, i.e. ``coeffs[0]`` is the
        constant term, ``coeffs[1]`` is the linear term, etc.

        Parameters
        ----------
        coeffs : iterable
            A tuple or list of `~astropy.units.Quantity` objects that have the
            appropriate velocity per time units.
        t0 : `astropy.time.Time`, numeric (optional)
            The reference time for the trend.
        """
        if coeffs is None:
            coeffs = []

        _coeffs = []

        for i, coeff in enumerate(coeffs):
            if not hasattr(coeff, 'unit'):
                raise ValueError("Input coefficients must be a Quantity with "
                                 "velocity / time^(-index) units!")

            if i == 0:
                self._v_unit = coeff.unit # record the unit of the constant term
            _unit = self._v_unit / (u.day ** i)

            if not coeff.unit.is_equivalent(_unit):
                raise u.UnitsError("Input coefficients must have velocity / "
                                   "time^(-index) units!")

            _coeffs.append(coeff.to(_unit))

        self.coeffs = tuple(_coeffs)
        self.t0 = t0

    def __call__(self, t):
        """Evaluate the predicted velocity at the specified time."""

        if len(self.coeffs) == 0:
            return np.zeros_like(t).astype(float)

        if isinstance(self.t0, Time):
            if not isinstance(t, Time):
                raise TypeError(
                    "t0 specified as an Astropy Time object, so you must pass "
                    "in an Astropy Time object here."
                )

        if self.t0 is not None:
            t = t - self.t0
            if hasattr(t, 'jd'):
                t = t.jd

        elif isinstance(t, Time):
            t = t.tcb.mjd

        coeffs = []
        for i, x in enumerate(self.coeffs):
            coeffs.append(x.to_value(self._v_unit / u.day**i))

        return np.polyval(coeffs[::-1], t) * self._v_unit

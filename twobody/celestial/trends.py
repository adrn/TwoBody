# Third-party
import astropy.units as u
import numpy as np

__all__ = ['PolynomialVelocityTrend']

class VelocityTrend(object):
    pass

class PolynomialVelocityTrend(VelocityTrend):
    """
    Represents a long-term radial velocity trend.

    Parameters
    ----------
    *coeffs :
        The coefficients of the polynomial in increasing order. That is, the
        power is the index of the coefficient, so the 0th coefficient is the
        constant, the 1st coefficient is the linear term, etc.

    TODO
    ----
    - We could add support for a reference epoch (other than mjd=0) in here...

    """
    def __init__(self, *coeffs):
        self.coeffs = []

        for i, coeff in enumerate(coeffs):
            if not hasattr(coeff, 'unit'):
                raise ValueError("Input coefficients must be a Quantity with "
                                 "velocity per time^(index) units!")

            if i == 0:
                self._v_unit = coeff.unit
            _unit = self._v_unit / (u.day ** i)

            if not coeff.unit.is_equivalent(_unit):
                raise u.UnitsError("Input coefficients must have velocity per "
                                   "time^(index) units!")

            self.coeffs.append(coeff.to(_unit))

        self.coeffs = list(coeffs)

    def __call__(self, t):
        if len(self.coeffs) == 0:
            raise ValueError("To evaluate the trend, you must have supplied "
                             "at least one coefficient value at creation.")
        t = np.atleast_1d(t)

        if not hasattr(t, 'unit'): # assume bare array has units = day
            t = t * u.day

        if t.unit.physical_type != 'time':
            raise u.UnitsError("Input time(s) must be a Quantity with time "
                               "units!")

        return np.polyval([x.value for x in self.coeffs[::-1]],
                          t.to(u.day).value) * self._v_unit

# class PiecewisePolynomialVelocityTrend(VelocityTrend):
#     # TODO:
#     # This class can also represent different, independent sections of the
#     # data to handle, e.g., calibration offsets between epochs. See the
#     # ``data_mask`` argument documentation below for more info.
#     pass

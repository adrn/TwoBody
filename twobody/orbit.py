# Third-party
import astropy.coordinates as coord
from astropy.coordinates.matrix_utilities import matrix_product
import astropy.units as u
import numpy as np
from numpy import pi

# Project
from . import elements as elem
from .transforms import get_t0, a1_sini, mf
from .anomaly import (eccentric_anomaly_from_mean_anomaly,
                      true_anomaly_from_eccentric_anomaly)
from .utils import _parse_time, ArrayProcessor
from .barycenter import Barycenter
from .wrap import cy_rv_from_elements


__all__ = ['KeplerOrbit']


class KeplerOrbit(object):

    def __init__(self, elements=None, elements_type='kepler',
                 barycenter=None, **kwargs):
        """

        Parameters
        ----------

        Examples
        --------


        """

        if elements is None:
            elements_cls = getattr(elem,
                                   "{0}Elements"
                                   .format(elements_type.capitalize()))

            # pass everything in kwargs to the class initializer
            elements = elements_cls(**kwargs)

        elif not isinstance(elements, elem.OrbitalElements):
            raise TypeError("'elements' must be an instance of an "
                            "OrbitalElements subclass.")

        if barycenter is not None and not isinstance(barycenter, Barycenter):
            raise TypeError("barycenter must be a twobody.Barycenter instance.")

        self.elements = elements
        self.barycenter = barycenter

    def __getattr__(self, name):

        if hasattr(self.elements, name):
            return getattr(self.elements, name)

        else:
            raise AttributeError("type object '{0}' has no attribute '{1}'"
                                 .format(self.__class__.__name__, name))

    def unscaled_radial_velocity(self, time,
                                 anomaly_tol=None, anomaly_maxiter=None):
        """Compute the radial velocity of the body at the specified times.

        Parameters
        ----------
        time : array_like, `astropy.time.Time`
            Array of times as barycentric MJD values, or an Astropy
            `~astropy.time.Time` object containing the times to evaluate at.
        anomaly_tol : numeric (optional)
            Tolerance passed to
            `~twobody.celestial.eccentric_anomaly_from_mean_anomaly` for solving
            for the eccentric anomaly. See default value in that function.
        anomaly_maxiter : numeric (optional)
            Maximum number of iterations to use in
            `~twobody.celestial.eccentric_anomaly_from_mean_anomaly` for solving
            for the eccentric anomaly. See default value in that function.

        Returns
        -------
        rv : numeric [m/s]
            Relative radial velocity - does not include systemtic velocity!
        """
        if anomaly_tol is None:
            # TODO: config item?
            anomaly_tol = 1E-10

        if anomaly_maxiter is None:
            # TODO: config item?
            anomaly_maxiter = 128

        # TODO: do we always want to use MJD? precision issues...
        time = time.tcb.mjd
        proc = ArrayProcessor(time)
        t, = proc.prepare_arrays()
        rv = cy_rv_from_elements(t, self.P.to(u.day).value, 1., self.e,
                                 self.omega.to(u.radian).value,
                                 self.M0.to(u.radian).value,
                                 self.t0.tcb.mjd,
                                 anomaly_tol, anomaly_maxiter)
        return np.atleast_1d(proc.prepare_result(rv))

    def radial_velocity(self, time):
        """Compute the ... TODO

        Parameters
        ----------
        time : array_like, `astropy.time.Time`
            Array of times as barycentric MJD values, or an Astropy
            `~astropy.time.Time` object containing the times to evaluate at.
        """

        # TODO: barycenter motion
        return self.K * self.unscaled_radial_velocity(time) # + self.barycenter.

    def reference_plane(self, time):
        """Compute the orbit at specified times in the two body barycentric
        frame aligned with the reference plane coordinate system (XYZ).

        Parameters
        ----------
        time : array_like, `astropy.time.Time`
            Array of times as barycentric MJD values, or an Astropy
            `~astropy.time.Time` object containing the times to evaluate at.
        """

        # mean anomaly
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            dt = time - self.t0
            M = 2*pi * dt / self.P - self.M0
            # TODO: should I do this:
            # M = coord.Angle(M.to(u.radian)).wrap_at(360*u.deg)
            M = M.to(u.radian)

        # eccentric anomaly
        E = eccentric_anomaly_from_mean_anomaly(M, self.e) * u.rad

        # true anomaly
        f = true_anomaly_from_eccentric_anomaly(E, self.e) * u.rad

        # distance from center of mass to orbiting body
        r = self.a * (1. - self.e * np.cos(E))

        # compute the orbit in the cartesian, orbital plane system (xyz):
        x = r * np.cos(f)
        y = r * np.sin(f)
        z = np.zeros_like(x)

        fac = 2*pi * self.a / self.P / np.sqrt(1 - self.e**2)
        vx = -fac * np.sin(f)
        vy = fac * (np.cos(f) + self.e)
        vz = np.zeros_like(vx)

        xyz = coord.CartesianRepresentation(x=x, y=y, z=z)
        vxyz = coord.CartesianDifferential(d_x=vx, d_y=vy, d_z=vz)

        # Construct rotation matrix to take the orbit from the orbital plane
        # system (xyz) to the reference plane system (XYZ):
        R1 = coord.matrix_utilities.rotation_matrix(self.omega, axis='z')
        R2 = coord.matrix_utilities.rotation_matrix(self.i, axis='x')
        R3 = coord.matrix_utilities.rotation_matrix(self.Omega, axis='z')
        Rot = R3 @ R2 @ R1

        # Rotate to the reference plane system
        XYZ = coord.CartesianRepresentation(matrix_product(Rot, xyz.xyz))
        VXYZ = coord.CartesianDifferential(matrix_product(Rot, vxyz.d_xyz))

        return XYZ.with_differentials(VXYZ)

    def icrs(self, time):
        pass

    def at_times(self, time):


        # Account for motion of the barycenter
        dt = times - self.barycenter.t0

        # TODO: transform
        frame = self.barycenter.reference_plane_frame
        bc = self.barycenter.frame.transform_to(reference_plane_frame)
        v_bary = bc.cartesian.differentials['s']

        XYZ = XYZ + bc.distance + v_bary * dt
        VXYZ = coord.CartesianDifferential(VXYZ.d_xyz + v_bary.d_xyz)


    # def t0(self, ref_mjd):
    #     """Un-mod the phase at pericenter, ``phi0``, to a time closest to the
    #     specified reference epoch.

    #     Parameters
    #     ----------
    #     ref_mjd : numeric
    #         Reference time in Barycentric MJD to get the pericenter time
    #         ``t0`` relative to.

    #     Returns
    #     -------
    #     t0 : `~astropy.time.Time`
    #         Pericenter time closest to input epoch.

    #     """
    #     return get_t0(self.phi0, self.P, ref_mjd)

    # def _generate_rv_curve(self, t):
    #     if isinstance(t, at.Time):
    #         _t = t.tcb.mjd
    #     else:
    #         _t = t

    #     rv = rv_from_elements(times=_t,
    #                           P=self.P.to(u.day).value,
    #                           K=self.K.to(u.m/u.s).value,
    #                           e=self.ecc,
    #                           omega=self.omega.to(u.radian).value,
    #                           phi0=self.phi0.to(u.radian).value,
    #                           anomaly_tol=self.anomaly_tol)

    #     if self.trend is not None:
    #         v_trend = self.trend(_t).to(u.m/u.s).value
    #         rv += v_trend

    #     return rv

    # def generate_rv_curve(self, t):
    #     """Generate a radial velocity curve evaluated at the specified times
    #     with the instantiated parameters.

    #     Parameters
    #     ----------
    #     t : array_like, `~astropy.time.Time`
    #         Time array. Either in BMJD or as an Astropy time.

    #     Returns
    #     -------
    #     rv : astropy.units.Quantity [speed]
    #     """
    #     rv = self._generate_rv_curve(t)
    #     return (rv*u.m/u.s).to(u.km/u.s)

    # def __call__(self, t):
    #     return self.generate_rv_curve(t)

    # def plot(self, t, ax=None, rv_unit=None, t_kwargs=None, plot_kwargs=None):
    #     """Plot the RV curve at the specified times.

    #     Parameters
    #     ----------
    #     t : array_like, `~astropy.time.Time`
    #         Time array. Either in BMJD or as an Astropy time.
    #     ax : `~matplotlib.axes.Axes`, optional
    #         The axis to draw on (default is to grab the current
    #         axes using `~matplotlib.pyplot.gca`).
    #     rv_unit : `~astropy.units.UnitBase`, optional
    #         Units to plot the radial velocities in (default is km/s).
    #     t_kwargs : dict, optional
    #         Keyword arguments passed to :class:`astropy.time.Time` with the
    #         input time array. For example, ``dict(format='mjd', scale='tcb')``
    #         for Barycentric MJD.
    #     plot_kwargs : dict, optional
    #         Any additional arguments or style settings passed to
    #         :func:`matplotlib.pyplot.plot`.

    #     Returns
    #     -------
    #     ax : `~matplotlib.axes.Axes`
    #         The matplotlib axes object that the RV curve was drawn on.

    #     """

    #     if ax is None:
    #         import matplotlib.pyplot as plt
    #         ax = plt.gca()

    #     if rv_unit is None:
    #         rv_unit = u.km/u.s

    #     if t_kwargs is None:
    #         t_kwargs = dict(format='mjd', scale='tcb')

    #     if plot_kwargs is None:
    #         plot_kwargs = dict()

    #     style = plot_kwargs.copy()
    #     style.setdefault('linestyle', '-')
    #     style.setdefault('alpha', 0.5)
    #     style.setdefault('marker', None)

    #     if not isinstance(t, at.Time):
    #         t = at.Time(t, **t_kwargs)
    #     rv = self.generate_rv_curve(t).to(rv_unit).value

    #     _t = getattr(getattr(t, t_kwargs['scale']), t_kwargs['format'])
    #     ax.plot(_t, rv, **style)

    #     return ax

    # # --------------------------------------------------------------------------
    # # Computed attributes

    # @property
    # def a1_sini(self):
    #     return a1_sini(self.P, self.K, self.ecc)

    # @property
    # def mf(self):
    #     return mf(self.P, self.K, self.ecc)

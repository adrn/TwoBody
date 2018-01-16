*********************
twobody Documentation
*********************

``twobody`` is a package for computing orbits and astronomical observables for
binary stars, exoplanets, and other gravitational two-body systems.

.. warning::

    TwoBody is still in beta. Use at your own risk!

For information about how to install the package see :ref:`install`. To learn
how and where to contribute, see the :ref:`contribute` page. For a quick
overview of available functionality, see the :ref:`getting-started` page. For a
quick primer on celestial mechanics and a description of the assumptions and
coordinate systems used in this package, see :ref:`celestial`. For more
exhaustive documentation on how to use the functions and classes in this
package, see the :ref:`full-api`.

.. toctree::
    :maxdepth: 1

    install
    contribute
    getting-started
    celestial
    full-api

~~~~~

.. plot::
    :context: close-figs
    :align: center
    :width: 512px

    from astropy.time import Time
    import astropy.units as u
    import astropy.coordinates as coord
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from twobody import KeplerOrbit, Barycenter

    origin = coord.ICRS(ra=14.745*u.deg, dec=71.512*u.deg,
                        distance=71.634*u.pc,
                        pm_ra_cosdec=32.123*u.mas/u.yr,
                        pm_dec=86.63*u.mas/u.yr,
                        radial_velocity=17.4123*u.km/u.s)
    b = Barycenter(origin=origin, t0=Time('J2000'))

    orb = KeplerOrbit(P=0.13*u.year, a=0.35*u.au, e=0.762,
                      omega=2.124*u.rad, M0=0*u.rad, t0=Time('J2000'),
                      Omega=121.53*u.deg, i=61*u.deg,
                      barycenter=b)

    t1 = Time('J2000') + np.linspace(0, orb.P.to(u.day).value, 10000)*u.day
    xyz = orb.orbital_plane(t1)
    XYZ = orb.reference_plane(t1).cartesian

    orbit_style = dict(marker='', linestyle='-', linewidth=2, color='#888888')
    body_style = dict(marker='o', linestyle='none', color='tab:red',
                      markersize=12, zorder=100)
    barycen_style = dict(marker='+', color='#888888', mew=2, ms=8)

    fig, axes = plt.subplots(1, 2, figsize=(8,4.5),
                             sharex=True, sharey=True)

    axes[0].plot(xyz.x, xyz.y, **orbit_style)
    axes[0].plot(0, 0, **barycen_style)

    axes[1].plot(XYZ.x, XYZ.y, **orbit_style)
    axes[1].plot(0, 0, **barycen_style)

    lim = 0.7
    for j in [0, 1]:
        axes[j].set_xlim(-lim, lim)
        axes[j].set_ylim(-lim, lim)


    axes[0].set_xlabel('$x$ [{0:latex_inline}]'.format(xyz.x.unit))
    axes[0].set_ylabel('$y$ [{0:latex_inline}]'.format(xyz.x.unit))

    axes[1].set_xlabel('$X$ [{0:latex_inline}]'.format(xyz.x.unit))
    axes[1].set_ylabel('$Y$ [{0:latex_inline}]'.format(xyz.x.unit))

    axes[0].set_title('orbital plane')
    axes[1].set_title('reference plane')

    fig.tight_layout()

    # --------------------
    # plot 2

    t2 = Time('J2000') + np.linspace(0, 3*orb.P.value, 10000)*orb.P.unit
    icrs = orb.icrs(t2)
    # TODO: this is only necessary because the released version of Astropy
    # doesn't support velocity transforms in SkyOffsetFrame
    _icrs = coord.ICRS(icrs.spherical.without_differentials())
    _origin = coord.ICRS(origin.spherical.without_differentials())
    offset = _icrs.transform_to(coord.SkyOffsetFrame(origin=_origin))
    # offset = icrs.transform_to(coord.SkyOffsetFrame(origin=origin))

    style = dict(marker='o', s=2, cmap='viridis')

    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5))

    axes[0].scatter(t2.mjd, icrs.radial_velocity,
                    c=t2.mjd, **style)
    axes[0].set_xlabel(r'time [${\rm MJD}$]')
    axes[0].set_ylabel(r'radial velocity [{0:latex_inline}]'
                       .format(icrs.radial_velocity.unit))
    axes[0].set_title('line of sight')

    axes[1].scatter(offset.lon.wrap_at(180*u.deg).milliarcsecond,
                    offset.lat.milliarcsecond,
                    c=t2.mjd, **style)
    axes[1].set_xlabel(r'$\Delta\alpha$ [{0:latex_inline}]'.format(u.mas))
    axes[1].set_ylabel(r'$\Delta\delta$ [{0:latex_inline}]'.format(u.mas))
    axes[1].set_title('sky offset')
    axes[1].set_xlim(-2, 37)
    axes[1].set_ylim(-2, 37)

    fig.tight_layout()


For developers
--------------

.. toctree::
    :maxdepth: 2

    tests

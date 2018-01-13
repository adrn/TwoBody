*********************
twobody Documentation
*********************

``twobody`` is a package for computing astronomical observables ...

.. toctree::
    :maxdepth: 1

    install
    contribute


Getting started
---------------

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from twobody import KeplerOrbit
    >>> import matplotlib.pyplot as plt
    >>> orb = KeplerOrbit(P=1.5*u.year, a=1.8*u.au, e=0.67,
    ...                   omega=17.14*u.deg, i=65*u.deg, Omega=5*u.deg,
    ...                   M0=0*u.deg, t0=Time('J2015.0'))
    >>> t = Time('2009-01-10') + np.linspace(0, 5, 1024) * u.year
    >>> rv = orb.radial_velocity(t)
    >>> plt.plot(t.datetime, rv.to(u.km/u.s), marker='') # doctest: +SKIP
    >>> plt.xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    >>> plt.ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP

.. plot::

    from astropy.time import Time
    import astropy.units as u
    from twobody import KeplerOrbit
    import matplotlib.pyplot as plt
    orb = KeplerOrbit(P=1.5*u.year, a=1.8*u.au, e=0.67,
                      omega=17.14*u.deg, i=65*u.deg, Omega=5*u.deg,
                      M0=0*u.deg, t0=Time('J2015.0'))
    t = Time('2009-01-10') + np.linspace(0, 5, 1024) * u.year
    rv = orb.radial_velocity(t)
    plt.plot(t.datetime, rv.to(u.km/u.s), marker='')
    plt.xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    plt.ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP


TODO: some examples

* Plot RV curve given elements
* Plot astrometric orbit given elements and barycenter motion

.. toctree::
    :maxdepth: 1

    celestial
    full-api


For developers
--------------

.. toctree::
    :maxdepth: 2

    tests

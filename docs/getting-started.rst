.. _getting-started:

********************************
Getting started with ``twobody``
********************************

For the examples below, we'll need the following imports:

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np


.. _getting-started-rv:

Compute and plot a radial velocity curve given orbital elements
---------------------------------------------------------------

We work with Keplerian orbits using the `~twobody.KeplerOrbit` class, which
accepts orbital elements on creation:

    >>> from twobody import KeplerOrbit
    >>> orb = KeplerOrbit(P=1.5*u.year, e=0.67,
    ...                   omega=17.14*u.deg, i=65*u.deg, Omega=0*u.deg,
    ...                   M0=35.824*u.deg, t0=Time('J2015.0'))

The elements are described in more detail in :ref:`celestial`, but, briefly,
``P`` is the orbital period, ``e`` is the eccentricity, ``omega`` is the
argument of pericenter, ``i`` is the inclination of the orbital plane, ``Omega``
is the longitude of the ascending node (which has no effect on the radial
velocity curve, so we set to 0), ``M0`` is the phase of the orbit at the
reference time ``t0``. All of the specific values of the above parameters were
made up at random. Note that we haven't specified a semi-major axis, so the
absolute scale of the radial velocity variations isn't specified. This is useful
if you are fitting for the velocity semi-amplitude ``K`` and systemic velocity
of the barycenter ``v0``, which are just linear parameters in the unscaled
radial velocity. We can compute the unscaled radial velocity by specifying an
array of times to the `~twobody.KeplerOrbit.unscaled_radial_velocity` method:

    >>> t = Time('2009-01-10') + np.linspace(0, 5, 1024) * u.year
    >>> unscaled_rv = orb.unscaled_radial_velocity(t)

Let's plot these values to visualize:

    >>> fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    >>> ax.plot(t.datetime, unscaled_rv.value, marker='') # doctest: +SKIP
    >>> ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    >>> ax.set_ylabel('unscaled radial velocity') # doctest: +SKIP

.. plot::
    :context: close-figs
    :align: center
    :width: 512px

    from astropy.time import Time
    import astropy.units as u
    import matplotlib.pyplot as plt
    import numpy as np
    from twobody import KeplerOrbit

    orb = KeplerOrbit(P=1.5*u.year, e=0.67,
                      omega=17.14*u.deg, i=65*u.deg, Omega=0*u.deg,
                      M0=35.824*u.deg, t0=Time('J2015.0'))

    t = Time('2009-01-10') + np.linspace(0, 5, 1024) * u.year
    unscaled_rv = orb.unscaled_radial_velocity(t)

    fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    ax.plot(t.datetime, unscaled_rv, marker='') # doctest: +SKIP
    ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    ax.set_ylabel('unscaled radial velocity') # doctest: +SKIP

If we do specify the semi-major axis, the velocity amplitude is known and we can
compute the actual (scaled) radial velocity curve with the
`~twobody.KeplerOrbit.radial_velocity` method:

    >>> orb = KeplerOrbit(P=1.5*u.year, e=0.67, a=1.77*u.au,
    ...                   omega=17.14*u.deg, i=65*u.deg, Omega=0*u.deg,
    ...                   M0=35.824*u.deg, t0=Time('J2015.0'))
    >>> fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    >>> ax.plot(t.datetime, rv.to(u.km/u.s).value, marker='') # doctest: +SKIP
    >>> ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    >>> ax.set_ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP

.. plot::
    :context: close-figs
    :align: center
    :width: 512px

    orb = KeplerOrbit(P=1.5*u.year, e=0.67, a=1.77*u.au,
                      omega=17.14*u.deg, i=65*u.deg, Omega=0*u.deg,
                      M0=35.824*u.deg, t0=Time('J2015.0'))
    rv = orb.radial_velocity(t)

    fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    ax.plot(t.datetime, rv.to(u.km/u.s).value, marker='') # doctest: +SKIP
    ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    ax.set_ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP

The radial velocities computed above are relative to the barycenter or reference
point of the orbit. We can add the barycentric or systemic velocity of the
system to the output radial velocities to get the actual line of sight
velocities. For example, if the systemic velocity is 100 km/s:

    >>> v0 = 100 * u.km/u.s
    >>> fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    >>> ax.plot(t.datetime, (rv + v0).to(u.km/u.s).value, marker='') # doctest: +SKIP
    >>> ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    >>> ax.set_ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP

.. plot::
    :context: close-figs
    :align: center
    :width: 512px

    v0 = 100 * u.km/u.s

    fig,ax = plt.subplots(1, 1) # doctest: +SKIP
    ax.plot(t.datetime, (rv + v0).to(u.km/u.s).value, marker='') # doctest: +SKIP
    ax.set_xlabel('time [{0:latex_inline}]'.format(u.year)) # doctest: +SKIP
    ax.set_ylabel('radial velocity [{0:latex_inline}]'.format(u.km/u.s)) # doctest: +SKIP

Note that both `~twobody.KeplerOrbit.radial_velocity` and
`~twobody.KeplerOrbit.unscaled_radial_velocity` assume that the barycenter does
not move tangentially when computing the velocity. That is, the radial velocity
computed is actually just :math:`\dot{Z}` computed in reference plane
coordinates (see :ref:`celestial`). For sources that move appreciable over the
baseline of observations, the observed line-of-sight velocity will change
slightly because of spherical projection effects, but the differences will be
small. See the docstring of `~twobody.KeplerOrbit.radial_velocity` for more
information.


.. _getting-started-astrometric:

Compute and plot an astrometric orbit curve given orbital elements and barycenter motion
----------------------------------------------------------------------------------------

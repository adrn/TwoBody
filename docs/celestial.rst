.. _celestial:

*****************************************************
Celestial mechanics and coordinate system conventions
*****************************************************

.. note::

    We mostly follow the notation and conventions defined in `Murray and Correia
    (2010) <https://arxiv.org/pdf/1009.1738.pdf>`_.

The Kepler ellipse
==================

We won't review the derivation or solution of the Kepler problem -- we recommend
sections 1–3 in `Murray and Correia (2010)
<https://arxiv.org/pdf/1009.1738.pdf>`_ -- and instead start from the fact that
the shape of a single Keplerian orbit is an ellipse. The position of an orbiting
body in its orbital plane is given by the vector :math:`\boldsymbol{r} =
\left(x, y, 0\right)`. The basis for this coordinate system is defined such that
:math:`\hat{x}` lies along the major axis of the orbit ellipse, increasing from
apocenter to pericenter. The basis vector :math:`\hat{z}` is aligned with the
orbital angular momentum and is perpendicular to :math:`\hat{x}`, and
:math:`\hat{y}` lies in the orbital plane and is perpendicular to both
:math:`\hat{x}` and :math:`\hat{z}`. Orbital motion is therefore constrained to
the :math:`x`-:math:`y` plane by construction. See figure below: the grey
ellipse represents the orbit of the (blue) body, with the arrow indicating the
direction of the orbit. The :math:`x`-:math:`y` basis is shown as dark black
arrows relative to the reference point (typically the center of mass):

.. plot::
    :width: 256px
    :align: center

    import numpy as np
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.pyplot as plt

    a = 1
    b = 0.75
    p = (b)**2 / a
    e = np.sqrt(1 - b**2/a**2)
    th = np.linspace(0, 2*np.pi, 1024)
    r = p / (1 + e*np.cos(th))

    x = r * np.cos(th)
    y = r * np.sin(th)

    fig, ax = plt.subplots(1, 1, figsize=(6,6))

    ax.plot(x, y, marker='', linestyle='-', linewidth=3, color='#aaaaaa')
    ax.plot(0, 0, marker='o')

    ax.plot(x[512], y[512], marker='o', markersize=10, color='tab:blue',
            markeredgewidth=1, markeredgecolor='k', zorder=100)
    ax.arrow(x[512], y[512], 0, -0.2, linewidth=3,
             zorder=10, head_width=0.04, head_length=0.05, color='tab:blue')

    ax.arrow(0, 0, 0, 1.5, color='k', zorder=-1, head_width=0.075, head_length=0.1)
    ax.arrow(0, 0, 1.5, 0, color='k', zorder=-1, head_width=0.075, head_length=0.1)
    ax.text(1.5, 0.1, r'$\hat{x}$', fontsize=24)
    ax.text(0.1, 1.5, r'$\hat{y}$', fontsize=24)

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

The position of the body on its orbit, :math:`\boldsymbol{r} = (x, y, 0)`, can
be expressed in terms of its orbital elements and time. First, we'll define the
angles typically used in expressing solutions to the Kepler problem. These are
the mean anomaly, :math:`M`, the eccentric anomaly, :math:`E`, and the true
anomaly, :math:`f`. These can be expressed in terms of the period of the orbit
:math:`P`, a time :math:`t`, a reference epoch :math:`t_0`, the orbital
eccentricity :math:`e`, and the semi-major axis :math:`a`:

.. math::

    M &= \frac{2\pi}{P} \, (t - t_0) - M_0 \\
    M &= E - e \, \sin{E} \\
    f &= 2 \, {\rm atan2}\left(\sqrt{1+e} \, \sin\frac{E}{2},
                               \sqrt{1-e} \, \cos\frac{E}{2}\right)\\
    r &= a \, (1 - e\,\cos{E})

In the above, :math:`r` is the distance of the body from the focus of the
ellipse closest to pericenter. :math:`M_0` is the mean anomaly or phase of the
orbit at the reference epoch. Using the above definitions, the position and
velocity of the body in the orbital plane coordinates is given by:

.. math::

    x &= r \, \cos{f} \\
    y &= r \, \sin{f} \\
    v_x &= \dot{r} \, \cos{f} - r \, \dot{f} \, \sin{f} \\
    &= -\frac{2\pi \, a}{P \, \sqrt{1 - e^2}} \, \sin{f} \\
    v_y &= \dot{r} \, \sin{f} + r \, \dot{f} \, \cos{f} \\
    &= \frac{2\pi \, a}{P \, \sqrt{1 - e^2}} \, \left[\cos{f} + e\right]

.. _celestial-reference-plane:

Observer or reference plane coordinates
=======================================

Of course, orbits of celestial bodies are generically rotated in with respect
to the observer's perspective (usually assumed to be sitting at the solar system
barycenter). The orientation of the orbit is therefore defined in terms of
orbital elements that rotate between the orbital plane :math:`(x, y, z)` system
and another reference coordinate system :math:`(X, Y, Z)`. In the new
:math:`(X, Y, Z)`, the reference plane must be defined. We use the tangent plane
at a point on the celestial sphere as the reference plane, with the observer
sitting along the positive :math:`Z` axis (see Figure 7 in `Murray and Correia
(2010) <https://arxiv.org/pdf/1009.1738.pdf>`_). At a given tangent point, the
:math:`\hat{X}` direction is aligned with North, and :math:`\hat{Y}` with the
East direction of the celestial coordinates used.

The angles that define the rotation from :math:`(x, y, z)` to :math:`(X, Y, Z)`
are the angular components of the orbital elements: longitude of the ascending
node, :math:`\Omega`, the argument of pericenter, :math:`\omega`, and the
inclination, :math:`i`. The full transformation is then a series of three
rotations: (1) rotate by :math:`\omega` around the :math:`z` axis to align
:math:`x'` with the line of nodes, (2) rotate by :math:`i` around :math:`x'`
to make the :math:`x'', y''` plane coincident with the reference plane, and (3)
rotate by :math:`\Omega` around :math:`z''` to align :math:`x'''` with
:math:`X`. So, to transform an orbit from its orbital plane to the reference
system, the full transformation is given by the composition of three rotation
matrices:

.. math::

    \begin{bmatrix} X \\ Y \\ Z \end{bmatrix} &=
        \boldsymbol{P}_{z}(\Omega) \,
        \boldsymbol{P}_{x}(i) \,
        \boldsymbol{P}_{z}(\omega) \,
        \begin{bmatrix} x \\ y \\ z \end{bmatrix}

where

.. math::

    \boldsymbol{P}_{x}(\phi) &=
        \begin{bmatrix}
            1 & 0 & 0 \\
            0 & \cos{\phi} & -\sin{\phi} \\
            0 & \sin{\phi} & \cos{\phi}
        \end{bmatrix} \\
    \boldsymbol{P}_{z}(\phi) &=
        \begin{bmatrix}
            \cos{\phi} & -\sin{\phi} & 0 \\
            \sin{\phi} & \cos{\phi} & 0 \\
            0 & 0 & 1
        \end{bmatrix}


See also:

* https://arxiv.org/pdf/1711.06601.pdf
* https://arxiv.org/pdf/1009.1738.pdf
* https://arxiv.org/pdf/1711.03595.pdf

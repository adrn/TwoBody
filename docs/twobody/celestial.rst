******************************************
Coordinate systems and celestial mechanics
******************************************

We mostly follow the notation and conventions defined in `Murray and Correia
(2010) <https://arxiv.org/pdf/1009.1738.pdf>`_. The position of an orbiting body
in its orbital plane is given by the vector :math:`\boldsymbol{r} = \left(x, y,
0\right)`. The basis for this coordinate system is defined such that
:math:`\hat{x}` lies along the major axis of the orbit ellipse, increasing from
apocenter to pericenter. The basis vector :math:`\hat{y}` also lies in the
orbital plane and is perpendicular to :math:`\hat{x}`. The :math:`\hat{z}`
direction is perpendicular to both, with the rotation of the full basis defined
such that (TODO: aligned with angular momentum of orbit, or???). Orbital motion
is therefore constrained to the :math:`x`-:math:`y` plane by construction.

Of course, orbits of celestial bodies are generically rotated in 3D with respect
to the observer's perspective (usually assumed to be sitting at the solar system
barycenter). The orbit is therefore typically defined in terms of orbital
elements that rotate between the :math:`(x, y, z)` system and some other
reference coordinate system :math:`(X, Y, Z)`. The reference plane must be
defined, but a typical use case in astronomy uses the tangent plane at a point
on the celestial sphere as the reference plane. The observer therefore sits
along the positive :math:`Z` axis. See Figure 7 in `Murray and Correia
(2010) <https://arxiv.org/pdf/1009.1738.pdf>`_. The angles that define the
rotation from :math:`(x, y, z)` to :math:`(X, Y, Z)` are the longitude of the
ascending node, :math:`\Omega`, the argument of pericenter, :math:`\omega`, and
the inclination, :math:`i`. TODO: use their definition of inclination, from 0 to
180 deg / prograde v. retrograde? The full transformation is then done in a
series of three rotations: (1) rotate by :math:`\omega` around the :math:`z`
axis to align :math:`x'` with the line of nodes, (2) rotate by :math:`i` around
:math:`x'` to make the :math:`x'', y''` plane coincident with the reference
plane, and (3) rotate by :math:`\Omega` around :math:`z''` to align :math:`x'''`
with :math:`X`. So, to transform an orbit from its orbital plane to the
reference system, the full transformation is given by the composition of three
rotation matrices:

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

TODO: transition...

The position of the body on its orbit, :math:`\boldsymbol{r} = (x, y, 0)`, can
be expressed in terms of orbital elements. First, some definitions:

.. math::

    M &= \frac{2\pi}{P} \, (t - t_0) \\
    M &= E - e \, \sin{E} \\
    f &= 2 \, {\rm atan2}\left(\sqrt{1+e} \, \sin\frac{E}{2},
                               \sqrt{1-e} \, \cos\frac{E}{2}\right)\\
    r &= a \, (1 - e\,\cos{E})

:math:`M` is the mean anomaly, :math:`t_0` is a reference time, :math:`E` is the
eccentric anomaly, :math:`f` is the true anomaly, and :math:`r` is the distance
of the body from the focus. Using the above definitions, the position and
velocity of the body in the orbital plane coordinates is:

.. math::

    x &= r \, \cos{f} \\
    y &= r \, \sin{f} \\
    v_x &= \dot{r} \, \cos{f} - r \, \dot{f} \, \sin{f} \\
    &= -\frac{2\pi \, a}{P \, \sqrt{1 - e^2}} \, \sin{f} \\
    v_y &= \dot{r} \, \sin{f} + r \, \dot{f} \, \cos{f} \\
    &= \frac{2\pi \, a}{P \, \sqrt{1 - e^2}} \, \left[\cos{f} + e\right]


See also:
* https://arxiv.org/pdf/1711.06601.pdf
* https://arxiv.org/pdf/1009.1738.pdf
* https://arxiv.org/pdf/1711.03595.pdf

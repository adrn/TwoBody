import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.time import Time
from twobody import KeplerOrbit, TwoBodyKeplerElements, ReferencePlaneFrame, Barycenter

elem = TwoBodyKeplerElements(P=1.5*u.year, e=0.67, M0=0*u.deg,
                             omega=17.14*u.deg, i=65*u.deg, Omega=0*u.deg,
                             m1=1*u.Msun, m2=2*u.Msun)

orb1 = KeplerOrbit(elem.primary)
orb2 = KeplerOrbit(elem.secondary)

t = Time('J2000') + np.linspace(0, 1, 10000) * elem.P

rp1 = orb1.reference_plane(t)
rp2 = orb2.reference_plane(t)

fig, ax = plt.subplots(1, 1, figsize=(5,5))

ax.plot(rp1.cartesian.x, rp1.cartesian.y,
        color='tab:blue', lw=10, marker='')
ax.plot(rp2.cartesian.x, rp2.cartesian.y,
        color='tab:red', lw=10, marker='')

i = 3100
ax.plot(rp1.cartesian.x[i], rp1.cartesian.y[i],
        marker='o', ms=75, color='tab:blue', mew=3, mec='w')
ax.plot(rp2.cartesian.x[i], rp2.cartesian.y[i],
        marker='o', ms=75, color='tab:red', mew=3, mec='w')

ax.scatter(0, 0, marker='o', s=80, color='k')

circ = mpl.patches.Circle((-0.5,0), 2.1, zorder=-100, color='w')
ax.add_patch(circ)

ax.set_aspect('equal', 'datalim')

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

ax.set_xlim(-2.6, 1.6)
ax.set_ylim(-2.1, 2.1)

for k in ax.spines:
    ax.spines[k].set_visible(False)

ax.set_facecolor('none')
fig.set_facecolor('none')

fig.tight_layout()

fig.savefig('icon.png', dpi=250, facecolor=fig.get_facecolor())
fig.savefig('icon.ico', format='png', dpi=50, facecolor=fig.get_facecolor())

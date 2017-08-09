"""
Celestial mechanics.

General comments
----------------
- Parameterization comes from Winn http://arxiv.org/abs/1001.2010
- Mean, eccentric, and true anomaly formulae from Wikipedia
  https://en.wikipedia.org/wiki/Eccentric_anomaly

"""

from .wrap import (cy_mean_anomaly_from_eccentric_anomaly,
                   cy_eccentric_anomaly_from_mean_anomaly_Newton1)

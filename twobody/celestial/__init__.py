"""
Celestial mechanics.

General comments
----------------
- Parameterization comes from Winn http://arxiv.org/abs/1001.2010
- Mean, eccentric, and true anomaly formulae from Wikipedia
  https://en.wikipedia.org/wiki/Eccentric_anomaly

"""

from .anomaly import *
from .core import *
from .orbit import *
from .trends import *

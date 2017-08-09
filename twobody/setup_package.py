from __future__ import absolute_import
from distutils.core import Extension
from astropy_helpers import setup_helpers

def get_extensions():
    exts = []

    # malloc
    # mac_incl_path = "/usr/include/malloc"

    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    # cfg['include_dirs'].append(mac_incl_path)
    cfg['include_dirs'].append('twobody/src')
    cfg['extra_compile_args'].append('--std=gnu99')
    cfg['sources'].append('twobody/src/twobody.c')
    exts.append(Extension('twobody.wrap', **cfg))

    return exts

def get_package_data():
    return {'twobody.wrap': ['src/*.h', 'src/*.c', '*.pyx']}

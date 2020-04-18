from collections import defaultdict
from distutils.core import Extension


def get_extensions():
    import numpy as np

    exts = []

    # malloc
    mac_incl_path = "/usr/include/malloc"

    cfg = defaultdict(list)
    cfg['include_dirs'].append(np.get_include())
    cfg['include_dirs'].append(mac_incl_path)
    cfg['include_dirs'].append('twobody/')
    cfg['extra_compile_args'].append('--std=gnu99')
    cfg['sources'].append('twobody/wrap.pyx')
    cfg['sources'].append('twobody/src/twobody.c')
    exts.append(Extension('twobody.wrap', **cfg))

    return exts


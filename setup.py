from collections import defaultdict

import numpy as np
from setuptools import Extension, setup

exts = []

cfg = defaultdict(list)
cfg["include_dirs"].append(np.get_include())
cfg["include_dirs"].append("twobody/")
cfg["extra_compile_args"].append("--std=gnu99")
cfg["sources"].append("twobody/wrap.pyx")
cfg["sources"].append("twobody/src/twobody.c")
exts.append(Extension("twobody.wrap", **cfg))

setup(ext_modules=exts)

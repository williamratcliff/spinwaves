from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("hello", ["hello.pyx", "McQueeny_Alg.c"])]

setup(
    name = 'McQueeny Cross Section Algorithm',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    include_dirs = [np.get_include()]
)
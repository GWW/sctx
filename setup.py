from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy, os
include_dirs = [numpy.get_include()]

extensions = [
    Extension("sctx.snvs.calc", ["src/sctx/snvs/calc.pyx"],
        include_dirs = include_dirs,
        language="c++",
        libraries=['pthread'],
        extra_compile_args=["-std=c++11", "-fPIC"],
        extra_link_args=["-std=c++11"],
    ),
]

pxd_dirs=['src/sctx/snvs/']
setup(name='sctx',
      packages=find_packages(where='src'),
      package_dir={"":"src"},

      ext_modules=cythonize(extensions, include_path=pxd_dirs),
)

from setuptools import setup
from Cython.Build import cythonize
import numpy, os, glob
from os.path import join
from setuptools.extension import Extension
include_dirs = [numpy.get_include()]

extensions = [
    Extension("sctx.snvs.calc", ["sctx/snvs/calc.pyx"],
        include_dirs = include_dirs,
        language="c++",
        libraries=['pthread'],
        extra_compile_args=["-std=c++11", "-fPIC"],
        extra_link_args=["-std=c++11"],
    ),
]

pxd_dirs=['sctx/snvs/']
setup(name='sctx',
      version='0.0',
      description='GW Single Cell Library',
      author='Gavin Wilson',
      author_email='gavin.w.wilson@gmail.com',
      license='MIT',
      packages=['sctx'],
      install_requires=['matplotlib','numpy','scipy', 'cython', 'sklearn', 
          'statsmodels', 'h5py', 'flammkuchen',
          'python-igraph', 'leidenalg', 'cyvcf2'],
      ext_modules=cythonize(extensions, include_path=pxd_dirs),
      scripts=['bin/sctxmisc'],
)

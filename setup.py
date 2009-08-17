from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
import os

module1 = Extension('rsa',
                    include_dirs = [],
                    libraries = ['gmp'],
                    library_dirs = [],
                    sources = ['src/python_rsa.cpp'])


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        module1,
    ]
)

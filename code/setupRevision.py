#!/usr/bin/env python
# encoding: utf-8

from numpy.distutils.core import setup, Extension

module1 = Extension('full_aep', sources=['full_aep.f90',
                                       'adStack.c',
                                       'adBuffer.f'],
                    extra_compile_args=['-O2', '-c'])

module2 = Extension('constraints', sources=['position_constraints.f90',
                                       'adStack.c',
                                       'adBuffer.f'],
                    extra_compile_args=['-O2', '-c'])

module3 = Extension('grid_param_fortran', sources=['grid_param.f90',
                                       'adStack.c',
                                       'adBuffer.f'],
                    extra_compile_args=['-O2', '-c'])

module4 = Extension('param_fortran', sources=['full_param.f90',
                                       'adStack.c',
                                       'adBuffer.f'],
                    extra_compile_args=['-O2', '-c'])

setup(
    version='0.0.1',
    description='parameterization paper code',
    ext_modules=[module1,module2,module3,module4],
    license='Apache License, Version 2.0',
)

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:13:44 2017

@author: Brittany
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
      ext_modules = cythonize("produceW.pyx"),
)
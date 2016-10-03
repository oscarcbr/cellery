#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  spearman.pyx part of cellery (ceRNAs linking inference)
#  
#  Copyright 2015 Oscar Bedoya Reina <oscarb@fgu1124.anat.ox.ac.uk>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

"""
Library to compile cspearman
"""

########################################################
#~ Import external libraries
########################################################
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np


########################################################
#~ Compile
########################################################
ext_modules = [Extension('cspearman',
    sources=["spearman.pyx", "cspearman.cc"],
    include_dirs = [np.get_include()],
    language="c++")]

setup(cmdclass         = {'build_ext': build_ext},
    ext_modules      = ext_modules)

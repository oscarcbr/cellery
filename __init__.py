#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  __init__.py part of cellery (ceRNAs linking inference)
#  
#  Copyright 2016 Oscar Bedoya Reina <obedoya@igmm-linux-005>
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

#NOTE: sqlite3 need to be compile with the option SQLITE_MAX_VARIABLE_NUMBER=3999

"""
Import all celery exceptions and packages
"""
__version__ = "1.01"

#----------------------------
#Import exceptions
from exceptions import CelleryStandardError,CelleryExceptionObjct, \
CelleryWarningObjct

#----------------------------
#Import external image libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#----------------------------
#Import external quantitative libraries
from itertools import product
import numpy as np
from numpy import apply_along_axis,divide,empty,flatnonzero,float64,inf, \
isfinite,linspace,log,ma,multiply,nan,random,warnings
import scipy as sp
from scipy import integrate,stats
from scipy.stats import norm,binom
from scipy.stats.mstats import normaltest,zmap

#----------------------------
#Import external text processing libraries
from string import lower,upper
from tempfile import NamedTemporaryFile
from Bio import Seq

#----------------------------
#Import external other os-wrapping libraries
import gzip
from multiprocessing import Queue,Process
import os
import shutil
from xml.dom import minidom
import BeautifulSoup

#----------------------------
#Import internal classes
from __class__ import gene
from __db__ import *

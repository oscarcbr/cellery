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
Test the correlation methods.
NOTE: To compile --> python setup.py build_ext -i
"""


########################################################
#~ Import libraries.
########################################################
cimport numpy as np
import cython
import numpy as np
import os
from numpy import float32,nan
from libc.math cimport isnan
from cpython cimport array
from libcpp.vector cimport vector
cdef extern double spearman_correlation(vector[double] a, vector[double] b)


########################################################
#~ Run Spearman's correlation 
########################################################
cpdef spearmanr_by(np.ndarray[float, ndim=1] a, \
	np.ndarray[float, ndim=1] b, int cntRplcts):
	"""
	Input: np.ndarray is the numpy array a. np.ndarray is the numpy 
	array b. c is the minimal number of time points/replicates to 
	calculate the correlation.
	Output: Spearman's correlation of a vs b.
	NOTE: np.nan and float values are valid.
	"""
	cdef lenA = len(a)
	cdef lenB = len(b)
	assert lenA==lenB
	cdef vector[double] aMskd
	cdef vector[double] bMskd
	cdef double aMskdPosA
	cdef double aMskdPosB
	cdef int rplcts = cntRplcts
	#look for non-nan values
	for posA in xrange(lenA):
		aMskdPosA = a[posA]
		aMskdPosB = b[posA]
		if not isnan(aMskdPosA) and not isnan(aMskdPosB):
			aMskd.push_back(aMskdPosA)
			bMskd.push_back(aMskdPosB)
			rplcts-=1
	#
	if rplcts>=0:
		return None
	else:
		return spearman_correlation(aMskd,bMskd)
		
########################################################
#~ Run Spearman's correlation for intervals
########################################################
def sprmnCrltnintrvls(aAVlsA,aAVlsB,intrvlsVlsA,intrvlsVlsB,
	tmpFldr,cntRplcts):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. intrvlsVlsA is an interval in aAVlsA 
	of 	interest to calculate Spearman's correlation. intrvlsVlsB is an 
	interval in aAVlsB of interest to calculate Spearman's correlation. 
	tmpFldr is a folder to store partial correlation results. cntRplcts 
	is the minimal number of time points/replicates to calculate the 
	correlation.
	Output: lSprmnVlsAVlsB is a list of variables A and B and Spearman's 
	correlations. If not present, the correlation file is written.
	"""
	intrvlA1,intrvlA2 = int(intrvlsVlsA[0]),int(intrvlsVlsA[1])
	intrvlB1,intrvlB2 = int(intrvlsVlsB[0]),int(intrvlsVlsB[1])
	outfl = os.path.join(tmpFldr,'%s_%s_%s_%s_%s.ecp'%(intrvlA1, \
	intrvlA2,intrvlB1,intrvlB2,cntRplcts))
	lSprmnVlsAVlsB = []
	if os.path.exists(outfl):
		for echl in open(outfl,'r'):
			echl = echl.strip().split()
			if echl:
				lSprmnVlsAVlsB.append((int(echl[0]),int(echl[1]), \
				float32(echl[2])))
	else:
		for vlsAPos in xrange(intrvlA1,intrvlA2):
			aVlsA = aAVlsA[vlsAPos]
			for vlsBPos in xrange(intrvlB1,intrvlB2):
				aVlsB = aAVlsB[vlsBPos]
				sprmnCrltn = spearmanr_by(aVlsA,aVlsB,cntRplcts)
				if sprmnCrltn != None:
					lSprmnVlsAVlsB.append((vlsAPos,vlsBPos, \
					float32(sprmnCrltn)))
		outfl = open(outfl,'w')
		outfl.write('\n'.join([' '.join([str(v) for v in vls]) for vls \
		in lSprmnVlsAVlsB]))
		outfl.close()
	return lSprmnVlsAVlsB

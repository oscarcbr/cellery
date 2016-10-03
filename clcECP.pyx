#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  clcECP.pyx part of cellery (ceRNAs linking inference)
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
Methods to calculated and manipulate correlations.
To compile:
cython -a clcECP.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /exports/igmm/software/pkg/el7/apps/python/2.7.10/include/python2.7 -o clcECP.so clcECP.c
"""


########################################################
#~ Import libraries.
########################################################
cimport cython
cimport numpy as np
import numpy as np
import os
from numpy import isnan as npisnan
from libc.math cimport isnan
from numpy cimport ndarray
from numpy import array,divide,float32,int8,ma,nan,zeros
from scipy import stats


########################################################
#~ Calculate ECP score for real values
########################################################
cpdef clcECP(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB):
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b.
	Output: Calculates ECP value (shared/total).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		vlObctA = aMrnVlsObctA[posMrn]
		vlObctB = aMrnVlsObctB[posMrn]
		total+=vlObctB
		total+=vlObctA
		if 	vlObctA>0 and vlObctB>0:
			shrd += vlObctB
			shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/total
		return ECP


########################################################
#~ Calculate ECP score for real values (including a masked miRNA array)
########################################################
cpdef clcECPMskd(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, np.ndarray \
	[np.int8_t, ndim=1] aMrnVlsMsk):
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. np.ndarray aMrnVlsMsk is a numpy array with bool
	values for masking aMrnVlsObctA and aMrnVlsObctB.
	Output: Calculates ECP value (shared/total).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		if aMrnVlsMsk[posMrn] == 0:#not masked
			vlObctA = aMrnVlsObctA[posMrn]
			vlObctB = aMrnVlsObctB[posMrn]
			total+=vlObctB
			total+=vlObctA
			if 	vlObctA>0 and vlObctB>0:
				shrd += vlObctB
				shrd += vlObctA
	#
	print shrd,total
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/total
		return ECP


########################################################
#~ Calculate ECP score density for real values
########################################################
cpdef clcECPDnsty(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, int lenA, int lenB):
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. int lenA is the length of ObctA. int lenB is the 
	length of ObctB.
	Output: Calculates ECP density value (shared/(total*(lenA+lenB))).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		vlObctA = aMrnVlsObctA[posMrn]
		vlObctB = aMrnVlsObctB[posMrn]
		total+=vlObctB
		total+=vlObctA
		if 	vlObctA>0 and vlObctB>0:
			shrd += vlObctB
			shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/(total*(lenA+lenB))
		return ECP


########################################################
#~ Calculate ECP score density for real values (including a masked miRNA
# array)
########################################################
cpdef clcECPDnstyMskd(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, np.ndarray \
	[np.int8_t, ndim=1] aMrnVlsMsk, int lenA, int lenB):
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. np.ndarray aMrnVlsMsk is a numpy array with bool 
	values for masking aMrnVlsObctA and aMrnVlsObctB. int lenA is the 
	length of ObctA.  int lenB is the length of ObctB.
	Output: Calculates ECP density value (shared/(total*(lenA+lenB))).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		if aMrnVlsMsk[posMrn] == 0:#not masked
			vlObctA = aMrnVlsObctA[posMrn]
			vlObctB = aMrnVlsObctB[posMrn]
			total+=vlObctB
			total+=vlObctA
			if 	vlObctA>0 and vlObctB>0:
				shrd += vlObctB
				shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/(total*(lenA+lenB))
		return ECP


########################################################
#~ Calculate ECP score for values excluding nan
########################################################
cpdef clcECPNan(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB):	
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b.
	Output: Calculates ECP value (sum(shared scores)/sum(all scores)).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		vlObctA = aMrnVlsObctA[posMrn]
		vlObctB = aMrnVlsObctB[posMrn]
		if not isnan(vlObctA):
			total += vlObctA
		if not isnan(vlObctB):
			total += vlObctB
		if 	not isnan(vlObctA) and not isnan(vlObctB):
			shrd += vlObctB
			shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/total
		return ECP


########################################################
#~ Calculate ECP score for values excluding nan (including a masked 
# miRNA array)
########################################################
cpdef clcECPNanMskd(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, np.ndarray \
	[np.int8_t, ndim=1] aMrnVlsMsk):	
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. np.ndarray aMrnVlsMsk is a numpy array with bool
	values for masking aMrnVlsObctA and aMrnVlsObctB.
	Output: Calculates ECP value (sum(shared scores)/sum(all scores)).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		if aMrnVlsMsk[posMrn] == 0:#not masked
			vlObctA = aMrnVlsObctA[posMrn]
			vlObctB = aMrnVlsObctB[posMrn]
			if not isnan(vlObctA):
				total += vlObctA
			if not isnan(vlObctB):
				total += vlObctB
			if 	not isnan(vlObctA) and not isnan(vlObctB):
				shrd += vlObctB
				shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/total
		return ECP


########################################################
#~ Calculate ECP score density for values excluding nan
########################################################
cpdef clcECPNanDnsty(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, int lenA, int lenB):	
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. int lenB is the length of ObctB. np.ndarray 
	aMrnVlsMsk is a numpy array with bool values for masking 
	aMrnVlsObctA and aMrnVlsObctB. int lenA is the length of ObctA. int 
	lenB is the length of ObctB.
	Output: Calculates ECP density value (shared/(total*(lenA+lenB))).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		vlObctA = aMrnVlsObctA[posMrn]
		vlObctB = aMrnVlsObctB[posMrn]
		if not isnan(vlObctA):
			total += vlObctA
		if not isnan(vlObctB):
			total += vlObctB
		if 	not isnan(vlObctA) and not isnan(vlObctB):
			shrd += vlObctB
			shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/(total*(lenA+lenB))
		return ECP


########################################################
#~ Calculate ECP score density for values excluding nan (and including a 
# masked miRNA array)
########################################################
cpdef clcECPNanDnstyMskd(np.ndarray[np.float32_t, ndim=1] aMrnVlsObctA, \
	np.ndarray[np.float32_t, ndim=1] aMrnVlsObctB, np.ndarray \
	[np.int8_t, ndim=1] aMrnVlsMsk, int lenA, int lenB):	
	"""
	Input: np.ndarray aMrnVlsObctA is the numpy array A. aMrnVlsObctB is 
	the numpy array b. int lenB is the length of ObctB. np.ndarray 
	aMrnVlsMsk is a numpy array with bool values for masking 
	aMrnVlsObctA and aMrnVlsObctB. int lenA is the length of ObctA. int 
	lenB is the length of ObctB.
	Output: Calculates ECP density value (shared/(total*(lenA+lenB))).
	"""
	cdef int posMrn = 0
	cdef double shrd = 0
	cdef double total = 0
	cdef double ECP = 0
	lenaMrnVlsObctB = len(aMrnVlsObctB)#assumes == len(aMrnVlsObctA)
	for posMrn in xrange(lenaMrnVlsObctB):
		if aMrnVlsMsk[posMrn] == 0:#not masked
			vlObctA = aMrnVlsObctA[posMrn]
			vlObctB = aMrnVlsObctB[posMrn]
			if not isnan(vlObctA):
				total += vlObctA
			if not isnan(vlObctB):
				total += vlObctB
			if 	not isnan(vlObctA) and not isnan(vlObctB):
				shrd += vlObctB
				shrd += vlObctA
	#
	if total == 0:
		return nan
	elif shrd == 0:
		return 0
	else:
		ECP = shrd/(total*(lenA+lenB))
		return ECP


########################################################
#~ Calculate ECP score for values excluding nan
########################################################
def rtrnECP(aAVlsA,aAVlsB,fldrOutECPPrws,intrvlsVlsA, \
	intrvlsVlsB,pntrCnts,aMskRef=False,aALenA=False,aALenB=False):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. intrvlsVlsA is an interval in aAVlsA 
	of interest to calculate ECP values. intrvlsVlsB is an interval in 
	aAVlsB of interest to calculate ECP values. fldrOutECPPrws is a 
	folder to store partial ECP value results. pntrCnts if True 
	indicates that aAVlsA and aAVlsB are counts so 0 values shall be 
	considered. aMskRef, aALenA, and aALenB are always False.
	Output: lECPVlsAVlsB is a list of variables A and B and ECP values 
	between them. If not present, the ECP value file is written.
	NOTE: For input arrays including nan values.
	"""
	#----------------------------
	#Find if input arrays has nan values
	if pntrCnts:
		ECPclcMthd = clcECP
	else:
		ECPclcMthd = clcECPNan
	#----------------------------
	#Run ECP calculation
	intrvlA1,intrvlA2 = int(intrvlsVlsA[0]),int(intrvlsVlsA[1])
	intrvlB1,intrvlB2 = int(intrvlsVlsB[0]),int(intrvlsVlsB[1])
	outfl = os.path.join(fldrOutECPPrws,'%s_%s_%s_%s.ecp'%(intrvlA1, \
	intrvlA2,intrvlB1,intrvlB2))
	lECPVlsAVlsB = []
	if os.path.exists(outfl):
		for echl in open(outfl,'r'):
			echl = echl.strip().split()
			if echl:
				lECPVlsAVlsB.append((int(echl[0]),int(echl[1]), \
				float32(echl[2])))
	else:
		for vlsAPos in xrange(intrvlA1,intrvlA2):
			aVlsA = aAVlsA[vlsAPos]
			for vlsBPos in xrange(intrvlB1,intrvlB2):
				aVlsB = aAVlsB[vlsBPos]
				ECPval = ECPclcMthd(aVlsA,aVlsB)
				if not isnan(ECPval):
					lECPVlsAVlsB.append((vlsAPos,vlsBPos, \
					float32(ECPval)))
		outfl = open(outfl,'w')
		outfl.write('\n'.join([' '.join([str(v) for v in vls]) for vls \
		in lECPVlsAVlsB]))
		outfl.close()
	return lECPVlsAVlsB


########################################################
#~ Calculate ECP score for values excluding nan (and including a masked 
# miRNA array)
########################################################			
def rtrnECPMskd(aAVlsA,aAVlsB,fldrOutECPPrws,intrvlsVlsA,intrvlsVlsB, \
	pntrCnts,aMskRef,aALenA=False,aALenB=False):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. intrvlsVlsA is an interval in aAVlsA 
	of interest to calculate ECP values. intrvlsVlsB is an interval in 
	aAVlsB of interest to calculate ECP values. fldrOutECPPrws is a 
	folder to store partial ECP value results. aMskRef is a mask array
	for the miRNAs (i.e. arrays within array A and B). pntrCnts if True 
	indicates that aAVlsA and aAVlsB are counts so 0 values shall be 
	considered. aALenA and aALenB are always False.
	Output: lECPVlsAVlsB is a list of variables A and B and ECP values 
	between them. If not present, the ECP value file is written.
	"""
	#----------------------------
	#Test input values and convert aMskRef for cython
	assert len(aAVlsA[0]) == aMskRef.size
	aMskRefCythn = aMskRef.view(dtype=int8)
	#----------------------------
	#Find if input arrays has nan values
	if pntrCnts:
		ECPclcMthd = clcECPMskd
	else:
		ECPclcMthd = clcECPNanMskd
	#----------------------------
	#Run ECP calculation
	intrvlA1,intrvlA2 = int(intrvlsVlsA[0]),int(intrvlsVlsA[1])
	intrvlB1,intrvlB2 = int(intrvlsVlsB[0]),int(intrvlsVlsB[1])
	outfl = os.path.join(fldrOutECPPrws,'%s_%s_%s_%s.ecp'%(intrvlA1, \
	intrvlA2,intrvlB1,intrvlB2))
	lECPVlsAVlsB = []
	if os.path.exists(outfl):
		for echl in open(outfl,'r'):
			echl = echl.strip().split()
			if echl:
				lECPVlsAVlsB.append((int(echl[0]),int(echl[1]), \
				float32(echl[2])))
	else:
		for vlsAPos in xrange(intrvlA1,intrvlA2):
			aVlsA = aAVlsA[vlsAPos]
			for vlsBPos in xrange(intrvlB1,intrvlB2):
				aVlsB = aAVlsB[vlsBPos]
				if not isnan(aMskRef[vlsAPos,vlsBPos]):
					ECPval = ECPclcMthd(aVlsA,aVlsB,aMskRefCythn)
					if not isnan(ECPval):
						lECPVlsAVlsB.append((vlsAPos,vlsBPos, \
						float32(ECPval)))
		outfl = open(outfl,'w')
		outfl.write('\n'.join([' '.join([str(v) for v in vls]) for vls \
		in lECPVlsAVlsB]))
		outfl.close()
	return lECPVlsAVlsB


########################################################
#~ Calculate ECP score density for values excluding nan
########################################################
def rtrnECPDnsty(aAVlsA,aAVlsB,fldrOutECPPrws,intrvlsVlsA, \
	intrvlsVlsB,pntrCnts,aMskRef=False,aALenA=False,aALenB=False):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. intrvlsVlsA is an interval in aAVlsA 
	of interest to calculate ECP values. intrvlsVlsB is an interval in 
	aAVlsB of interest to calculate ECP values. fldrOutECPPrws is a 
	folder to store partial ECP value results. pntrCnts if True 
	indicates that aAVlsA and aAVlsB are counts so 0 values shall be 
	considered. aALenA is an array of object lengths in the same order 
	that aAVlsA. aALenB is an array of object lengths in the same order 
	that aAVlsB. aMskRef is always False.
	Output: lECPVlsAVlsB is a list of variables A and B and ECP values 
	between them. If not present, the ECP value file is written.
	"""
	#----------------------------
	#Test input values
	assert aALenA and aALenB
	assert len(aALenA)==len(aAVlsA)
	assert len(aALenB)==len(aAVlsB)
	#----------------------------
	#Find if input arrays has nan values
	if pntrCnts:
		ECPclcMthd = clcECPDnsty
	else:
		ECPclcMthd = clcECPNanDnsty
	#----------------------------
	#Run ECP calculation
	intrvlA1,intrvlA2 = int(intrvlsVlsA[0]),int(intrvlsVlsA[1])
	intrvlB1,intrvlB2 = int(intrvlsVlsB[0]),int(intrvlsVlsB[1])
	outfl = os.path.join(fldrOutECPPrws,'%s_%s_%s_%s.ecp'%(intrvlA1, \
	intrvlA2,intrvlB1,intrvlB2))
	lECPVlsAVlsB = []
	if os.path.exists(outfl):
		for echl in open(outfl,'r'):
			echl = echl.strip().split()
			if echl:
				lECPVlsAVlsB.append((int(echl[0]),int(echl[1]), \
				float32(echl[2])))
	else:
		for vlsAPos in xrange(intrvlA1,intrvlA2):
			aVlsA = aAVlsA[vlsAPos]
			lenA = aALenA[vlsAPos]
			for vlsBPos in xrange(intrvlB1,intrvlB2):
				aVlsB = aAVlsB[vlsBPos]
				lenB = aALenB[vlsBPos]
				ECPval = ECPclcMthd(aVlsA,aVlsB,lenA,lenB)
				if not isnan(ECPval):
					lECPVlsAVlsB.append((vlsAPos,vlsBPos, \
					float32(ECPval)))
		outfl = open(outfl,'w')
		outfl.write('\n'.join([' '.join([str(v) for v in vls]) for vls \
		in lECPVlsAVlsB]))
		outfl.close()
	return lECPVlsAVlsB


########################################################
#~ Calculate ECP score density for values excluding nan for masked 
# arrays
########################################################			
def rtrnECPDnstyMskd(aAVlsA,aAVlsB,fldrOutECPPrws,intrvlsVlsA, \
	intrvlsVlsB,pntrCnts,aMskRef=False,aALenA=False,aALenB=False):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. intrvlsVlsA is an interval in aAVlsA 
	of interest to calculate ECP values. intrvlsVlsB is an interval in 
	aAVlsB of interest to calculate ECP values. fldrOutECPPrws is a 
	folder to store partial ECP value results. pntrCnts if True 
	indicates that aAVlsA and aAVlsB are counts so 0 values shall be 
	considered. aMskRef is a mask array for the miRNAs (i.e. arrays 
	within array A and B). aALenA is an array of object lengths in the 
	same order that aAVlsA. aALenB is an array of object lengths in the 
	same order that aAVlsB.
	Output: lECPVlsAVlsB is a list of variables A and B and ECP values 
	between them. If not present, the ECP value file is written.
	"""
	#----------------------------
	#Test input values
	assert len(aALenA)==len(aAVlsA)
	assert len(aALenB)==len(aAVlsB)
	assert len(aAVlsA[0]) == aMskRef.size
	aMskRefCythn = aMskRef.view(dtype=int8)
	#----------------------------
	#Find if input arrays has nan values
	if pntrCnts:
		ECPclcMthd = clcECPDnstyMskd
	else:
		ECPclcMthd = clcECPNanDnstyMskd
	#----------------------------
	#Run ECP calculation
	intrvlA1,intrvlA2 = int(intrvlsVlsA[0]),int(intrvlsVlsA[1])
	intrvlB1,intrvlB2 = int(intrvlsVlsB[0]),int(intrvlsVlsB[1])
	outfl = os.path.join(fldrOutECPPrws,'%s_%s_%s_%s.ecp'%(intrvlA1, \
	intrvlA2,intrvlB1,intrvlB2))
	lECPVlsAVlsB = []
	if os.path.exists(outfl):
		for echl in open(outfl,'r'):
			echl = echl.strip().split()
			if echl:
				lECPVlsAVlsB.append((int(echl[0]),int(echl[1]), \
				float32(echl[2])))
	else:
		for vlsAPos in xrange(intrvlA1,intrvlA2):
			aVlsA = aAVlsA[vlsAPos]
			lenA = aALenA[vlsAPos]
			for vlsBPos in xrange(intrvlB1,intrvlB2):
				aVlsB = aAVlsB[vlsBPos]
				lenB = aALenB[vlsBPos]
				if not isnan(aMskRef[vlsAPos,vlsBPos]):
					ECPval = ECPclcMthd(aVlsA,aVlsB,aMskRefCythn,lenA, \
					lenB)
					if not isnan(ECPval):
						lECPVlsAVlsB.append((vlsAPos,vlsBPos, \
						float32(ECPval)))
		outfl = open(outfl,'w')
		outfl.write('\n'.join([' '.join([str(v) for v in vls]) for vls \
		in lECPVlsAVlsB]))
		outfl.close()
	return lECPVlsAVlsB

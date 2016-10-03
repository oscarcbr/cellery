#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  matrices.py part of cellery (ceRNAs linking inference)
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

"""
Libraries to manipulate and merge matrices and arrays
"""

########################################################
#~ Import libraries.
########################################################
from copy import copy
from numpy import array,float32,float64,empty,ma,nan,vectorize

import numpy as np


########################################################
#~ Return a mask for isoforms coming from the same gene
########################################################
def mkMskIsfrmsSmGn(aVlsAVlsB,aVlsAIsfNms,aVlsBIsfNms,dIsfNmsGnNms):
	"""
	Input: aVlsAVlsB is a matrix of size valsA x valsB with values on 
	isoforms in set A and B respectively. aVlsAIsfNms is an array with
	isoform names in the same order as rows in aVlsAVlsB. aVlsBIsfNms 
	is an array with isoform names in the same order as columns in 
	aVlsAVlsB. dIsfNmsGnNms is a dictionary 
	Output: aMskIsfrmsSmGn is an array of size valsA x valsB  with 
	boolean values representing a masked value in case isoforms for 
	valA row and valB column are coming from the same gene.
	"""
	#--------------------------
	# test input variables
	sIsfmNms = set(dIsfNmsGnNms.keys())
	lenaVlsAIsfNms,lenaVlsBIsfNms = len(aVlsAIsfNms),len(aVlsBIsfNms)
	assert aVlsAVlsB.shape = (lenaVlsAIsfNms,lenaVlsBIsfNms)
	assert set(aVlsAIsfNms).issubset(sIsfmNms) and set(aVlsBIsfNms). \
	issubset(sIsfmNms)
	#--------------------------
	# make output mask matrix and fill it
	aMskIsfrmsSmGn = np.zeros(aVlsAVlsB.shape,dtype=bool)
	for isfrmAPos in xrange(lenaVlsAIsfNms):
		isfrmANm = aVlsAIsfNms[isfrmAPos]
		gnNmA = dIsfNmsGnNms[isfrmANm]
		for isfrmBPos in xrange(lenaVlsBIsfNms):
			isfrmBNm = aVlsBIsfNms[isfrmBPos]
			gnNmB = dIsfNmsGnNms[isfrmBNm]
			if gnNmA==gnNmB:
				aMskIsfrmsSmGn[isfrmAPos][isfrmBPos] = True
	#--------------------------
	# return mask
	return aMskIsfrmsSmGn


########################################################
#~ Return a mask for genes/isoforms coming from the same gene
########################################################
def mkMskNonOvrlpng(aVlsAVlsB,aVlsAChr,aVlsBChr,aVlsAMinInIntrvls, \
	aVlsBMinInIntrvls,aVlsAMaxInIntrvls,aVlsBMaxInIntrvls,N=10000):
	"""
	Input: aVlsAVlsB is a matrix of size valsA x valsB with values on 
	isoforms in set A and B respectively. aVlsAChr is an array with
	chromosomes in the same order as rows in aVlsAVlsB. aVlsBChr is an 
	array with chromosomes in the same order as columns in aVlsAVlsB. 
	aVlsAMinInIntrvls is an array with the minimum isoform/gene position 
	in the same order as rows in aVlsAVlsB. aVlsBMinInIntrvls is an 
	array with the minimum isoform/gene position in the same order as 
	columns in aVlsAVlsB. aVlsAMaxInIntrvls is an array with the maximum 
	isoform/gene position in the same order as rows in aVlsAVlsB. 
	aVlsBMaxInIntrvls is an array with the maximum isoform/gene position 
	in the same order as columns in aVlsAVlsB. N is the distance in bps
	under which isoforms or genes are going to be masked.
	Output: aMskNonOvrlpng is an array of size valsA x valsB  with 
	boolean values representing a masked value in case genes/isoforms 
	for valA row and valB column are in proximity by at most N.
	"""
	#--------------------------
	# Method to define overlaping
	def overlap(start1, end1, start2, end2):
		return end1 >= start2 and end2 >= start1
	#--------------------------
	# test input variables
	lenaVlsAChr,lenaVlsBChr = len(aVlsAChr),len(aVlsBChr)
	lenaVlsAMinInIntrvls,lenaVlsBMinInIntrvls = len(aVlsAMinInIntrvls), \
	len(aVlsBMinInIntrvls)
	lenaVlsAMaxInIntrvls,lenaVlsBMaxInIntrvls = len(aVlsAMaxInIntrvls), \
	len(aVlsBMaxInIntrvls)
	assert aVlsAVlsB.shape = (lenaVlsAChr,lenaVlsBChr)
	assert lenaVlsAChr==lenaVlsAMinInIntrvls and \
	lenaVlsBChr==lenaVlsBMinInIntrvls
	assert lenaVlsAChr==lenaVlsAMaxInIntrvls and \
	lenaVlsBChr==lenaVlsBMaxInIntrvls
	#--------------------------
	# Mask isoforms/genes in proximity for N/2 nts in each direction,
	# adding N in total.
	half = N/2#to add for each gene. and N must be a paired number
	aMskNonOvrlpng = np.zeros(aVlsAVlsB.shape,dtype=bool)
	for isfrmAPos in xrange(lenaVlsAChr):
		isfrmAChr = aVlsAChr[isfrmAPos]
		isfrmAMinPos = aVlsAMinInIntrvls[isfrmAPos]-half
		isfrmAMaxPos = aVlsAMaxInIntrvls[isfrmAPos]+half
		for isfrmBPos in xrange(lenaVlsBChr):
			isfrmBChr = aVlsBChr[isfrmBPos]
			isfrmBMinPos = aVlsBMinInIntrvls[isfrmBPos]-half
			isfrmBMaxPos = aVlsBMaxInIntrvls[isfrmBPos]+half
			if isfrmAChr==isfrmBChr and overlap(isfrmAMinPos, \
				isfrmAMaxPos,isfrmBMinPos,isfrmBMaxPos):
				aMskNonOvrlpng[isfrmAPos][isfrmBPos] = True
	#--------------------------
	# return mask	
	return aMskNonOvrlpng


########################################################
#~ Merge three large matrices (AxA, AxB, and BxB) into a single one.
########################################################
def mrgThreeLrgaMetrics(aMetrcAxA,aMetrcBxB,aMetrcAxB):
	"""
	Input: Three different matrices {aMetrcAxA,aMetrcBxB,aMetrcAxB}. 
	aMetrcAxA has dimensions AxA, aMetrcBxB has dimensions BxB, and 
	aMetrcAxB has dimension AxB.
	Output:  aJndMtrc is a joined metric with aMetrcAxA in the top left
	corner, aMetrcBxB in the bottom right corner, aMetrcAxB in the top
	right corner, and aMetrcAxB.T in the bottom left corner.
	"""
	aJndMtrc = np.append(aMetrcAxA,aMetrcAxB,axis=1)
	bttm = np.append(aMetrcAxB.T,aMetrcBxB,axis=1)
	aJndMtrc = np.append(aJndMtrc,bttm,axis=0)
	return aJndMtrc


########################################################
#~ Merge two objects into a single one. 
########################################################
def mrglObjcts(lObjctsA,lObjctsB):
	"""
	Input: Two list of objects to be merged into a single one {lObjctsA 
	and lObjctsB}. 
	Output: lMrgdObjct is a list of object of length len(lObjctsA)+
	len(lObjctsB), with objects in order lObjctsA and further lObjctsB.
	The objects in lMrgdObjct have updated position in order lObjctsA+
	lObjctsB (i.e. objct.pos in lObjctsB = objct.pos+len(lObjctsA)).
	"""
	#--------------------------
	# Copy all objects in input order and update their pos
	lMrgdObjct = []
	pos = -1
	for objct in lObjctsA:
		pos+=1
		cpyObjct = copy(objct)
		cpyObjct.pos = pos
		lMrgdObjct.append(cpyObjct)
	for objct in lObjctsB:
		pos+=1
		cpyObjct = copy(objct)
		cpyObjct.pos = pos
		lMrgdObjct.append(copy(objct))
	return lMrgdObjct


########################################################
#~ Function to be vectorized to calculate average values for matrix 
# (i.e. per gene) given a supermatrix and a set of positions
########################################################
def rtrnAvrgMskdArray(aValsAValsBMskd,aRowsPosValsA,aClmnsPosValsB):
	"""
	Input: aValsAValsBMskd is a masked array (the supermatrix) with 
	background metric values to be averaged following aRowsPosValsA and 
	aClmnsPosValsB. aRowsPosValsA is an array with the position of rows 
	in the supermatrix	(aValsAValsBMskd). aClmnsPosValsB is an array 
	with the positions of columns in the supermatrix.
	Output: avrgMskdVal is the average of the masked substracted 
	positions in aRowsPosValsA and aClmnsPosValsB of the supermatrix 
	(aValsAValsBMskd).
	"""
	if len(aRowsPosValsA) and len(aClmnsPosValsB):
		mskdECPsPerGn = ma.mean(aValsAValsBMskd[aRowsPosValsA,:] \
		[:,aClmnsPosValsB]).__float__()
		avrgMskdVal = float32(mskdECPsPerGn)
		return avrgMskdVal
	else:
		return nan


########################################################
#~ Return an array with positions of all isoform for each genes.
########################################################
def rtrnaObjctGnPos(lObjct,dObjctGn):
	"""
	Input: lObjct is a list of objects with pointer "name". dObjctGn is 
	a dictionary of isoforms to gene names.
	Output: aObjctOut is an array with sorted gene names and positions of
	isoforms in lObjct as sub-arrays. dSrtdObjctCdPos is the gene names 
	and their positions in aObjctOut. dIsfrmsPosGnsPos are the isoforms 
	in lObjct and their corresponding gene position in aObjctOut.
	NOTE: Returns an array of arrays in which each position is a gene 
	array with the positions of different isoforms within it.
	"""
	srtdObjctCd = sorted(set([dObjctGn[objct.name] for objct in \
	lObjct]))
	dSrtdObjctCdPos = dict([(gnNm,pos) for pos,gnNm in \
	enumerate(srtdObjctCd)])
	dIsfrmsPosGnsPos = {}
	#--------------------------
	#create an empty list of list to store the positions of the isoforms
	aObjctOut = [[] for x in xrange(len(srtdObjctCd))]
	for objct in lObjct:
		name = objct.name
		pos = objct.pos
		gnNm = dObjctGn[name]
		gnPos = dSrtdObjctCdPos[gnNm]
		aObjctOut[gnPos].append(pos)
		dIsfrmsPosGnsPos[pos] = gnPos
	aObjctOut = array(aObjctOut)
	return aObjctOut,dSrtdObjctCdPos,dIsfrmsPosGnsPos


########################################################
#~ Return the intersection of two matrices of the same size
########################################################
def rtrnIntrsctnTwoMtrx(aVlsAVlsB1,aVlsAVlsB2,vrbse=True, \
	ststc='multiply'):
	"""
	Input: Two different matrices with the same dimensions {aVlsAVlsB1,
	aVlsAVlsB2}. If vrbse is True, the log for the run is going 
	to be printed.
	Output:  aVlsAVlsBI is an intersecting matrix, where all nan values
	are going to be nan, and the intersecting cells will have a value of
	ststc
	"""
	aVlsAVlsB1Msk = ma.masked_invalid(aVlsAVlsB1)
	aVlsAVlsB2Msk = ma.masked_invalid(aVlsAVlsB2)
	aVlsAVlsBI = getattr(ma,ststc)(aVlsAVlsB1Msk,aVlsAVlsB2Msk)
	#----------------------------
	# report
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		(aVlsAVlsBI.count(),'edges were found to be shared')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	# fill output matrix
	aVlsAVlsBI.fill_value = nan
	aVlsAVlsBI = aVlsAVlsBI.filled()
	return aVlsAVlsBI,ovrAllogRun


########################################################
#~ Return a mask for rows and columns with less than "n" degree 
# values/nodes
########################################################
def rtrnMskDegreeGrtrN(aMsk,n):
	"""
	Input: aMsk is a binary mask matrix. n is the number of degrees of 
	interest (under which values are going to be masked).
	Output: aMskDegreeGrtrN is a matrix where all nodes that have a 
	degree equal or higher than n with 1 values and the other with 0.
	NOTE: Returns a masking array with 1 and 0 values. It masks and 
	replaces an input matrix of binary data with nodes with degrees less 
	than n replaced by 1 (i.e. masked). It runs iteratively until 
	converging (all nodes have a degree equal or higher than n).
	"""
	#--------------------------
	# invert mask
	aMskDegreeGrtrN = np.logical_not(aMsk)
	#--------------------------
	# calculate degree of nodes
	rowSum = np.sum(aMskDegreeGrtrN,1)
	clmnSum = np.sum(aMskDegreeGrtrN,0)
	crrntCntRow = len(rowSum)#to start no mask
	crrntCntClmn = len(clmnSum)#to start no mask
	#--------------------------
	# mask nodes with less degree than n
	mksdRowSum = ma.masked_less(rowSum,n)
	mksdClmnSum = ma.masked_less(clmnSum,n)
	cntMskdRow = mksdRowSum.count()
	cntMskdClmn = mksdClmnSum.count()
	if not cntMskdRow or not cntMskdClmn:#no answer
		return False
	#--------------------------
	# make degree mask to converge
	while (cntMskdRow!=crrntCntRow or cntMskdClmn!=crrntCntClmn) and \
		cntMskdRow>0 and cntMskdClmn>0:#test for covergence
		crrntCntRow = cntMskdRow
		crrntCntClmn = cntMskdClmn
		posToZeroRow = np.where(mksdRowSum.mask==1)[0]
		posToZeroClmn = np.where(mksdClmnSum.mask==1)[0]
		if len(posToZeroRow)>0:
			aMskDegreeGrtrN[posToZeroRow,:] = 0
		if len(posToZeroClmn)>0:
			aMskDegreeGrtrN[:,posToZeroClmn] = 0
		rowSum = np.sum(aMskDegreeGrtrN,1)
		clmnSum = np.sum(aMskDegreeGrtrN,0)
		mksdRowSum = ma.masked_less(rowSum,n)
		mksdClmnSum = ma.masked_less(clmnSum,n)
		cntMskdRow = mksdRowSum.count()
		cntMskdClmn = mksdClmnSum.count()
		if not cntMskdRow or not cntMskdClmn:#no answer
			return False
	#--------------------------
	# invert output mask so nodes with higher degrees than n be unmasked
	aMskDegreeGrtrN = np.logical_not(aMskDegreeGrtrN)
	return aMskDegreeGrtrN


########################################################
#~ Return a mask for rows and columns with less than "n" degree 
# values/nodes. The input are three matrices: ECP values, a 
# same-gene-isoform mask, and a valid expression mask.
########################################################
def unifyMsks3Mtrx(aMetrcECPflld,aMskIsfrmsSmGn,aMskdExpr,minDgr):
	"""
	Input: aMetrcECPflld is a matrix with ECP values or nan where 
	0|null. aMskIsfrmsSmGn is a matrix with a mask (0,1) for genes 
	coming from the same gene. aMskdExpr is a matrix with a mask for 
	genes whose correlation in expression is of interest (positive for 
	instance). minDgr is a degree number.
	Output: aMskDegreeGrtrN is the unified mask with only pairs 
	following the previous description.
	NOTE: Unify masks from three different sources. The first one is a 
	matrix with ECP values with nan in place of 0|null. The other two 
	matrices are masks with binary data. It retrieves a unified mask for 
	genes with degrees higher than n.
	"""
	aMsk = ma.masked_invalid(aMetrcECPflld).mask
	aMsk = ma.mask_or(aMsk,aMskIsfrmsSmGn)
	aMsk = ma.mask_or(aMsk,aMskdExpr)
	#--------------------------
	# mask node with nodes of a degree lower than minDgr
	aMskDegreeGrtrN = rtrnMskDegreeGrtrN(aMsk,minDgr)
	return aMskDegreeGrtrN


########################################################
#~ Wrapper to calculate average values for matrix (i.e. per gene) given 
# a supermatrix and a set of positions
########################################################
def wrprRtrnMskdAvrgECPGnLncrn(aValsAValsBMskd,aRowsPosValsA, \
	aClmnsPosValsB):
	"""
	Input: aValsAValsBMskd is a masked array (the supermatrix) with 
	background metric values to be averaged following aRowsPosValsA and 
	aClmnsPosValsB. aRowsPosValsA is an array with the position of rows 
	in the supermatrix	(aValsAValsBMskd). aClmnsPosValsB is an array 
	with the positions of columns in the supermatrix.
	Output: mskdValsAValsBPerGn is the average of the masked substracted 
	positions in aRowsPosValsA and aClmnsPosValsB of the supermatrix 
	(aValsAValsBMskd).
	NOTE: Is a wrapper to return the masked calculations using
	rtrnAvrgMskdArray.
	"""
	vctrzRtrnAvrgMskdArray = vectorize(rtrnAvrgMskdArray,excluded={0})
	mskdValsAValsBPerGn = vctrzRtrnAvrgMskdArray(aValsAValsBMskd, \
	aRowsPosValsA,aClmnsPosValsB)
	return mskdValsAValsBPerGn

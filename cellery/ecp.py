#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ecp.py part of cellery (ceRNAs linking inference)
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
Methods to calculate ECP values (Endogenous Competition Potential)
"""

########################################################
#~ Import libraries.
########################################################
from cellery import exceptions
from itertools import product
from multiprocessing import Queue,Process
from numpy import array,empty,float32,float64,nan,zeros
from clcECP import rtrnECP,rtrnECPMskd,rtrnECPDnsty,rtrnECPDnstyMskd

import os
import sqlite3


########################################################
#~ Compute ECP values for all combinations of two arrays of arrays with 
# values.
########################################################
def cmpECP(aMrnVlsDtA,aMrnVlsDtB,aANmsA,aANmsB,fldrOutECPPrws, \
	aALenA=False,aALenB=False,aMskRef=False,nThrds=10,intrvlSz=700, \
	sqlFl=False,pntrCnts=True):
	"""
	Input: aMrnVlsDtA is an array A of arrays with values for miRNAs.
	aMrnVlsDtB is an array B of arrays with values for miRNAs. aANmsA is 
	the array of variable names in the same position as the numbers in 
	vrblAPos. aANmsB is the array of variable names in the same order as 
	vrblBPos. fldrOutECPPrws is a folder to store partial ECP results. 
	Optionally, aALenA is an array of object lengths in the same order 
	that aAVlsA. aALenB is an array of object lengths in the same order 
	that aAVlsB. aMskRef is a mask array for the miRNAs (i.e. arrays 
	within array A and B). nThrds is the number of threads to run in 
	parallel. intrvlSz is the size of the interval to run in multithread. 
	sqlFl is a sql database to save the ECP values. If pntrCnts is True 
	aAVlsA and aAVlsB are counts so 0 values shall be considered 
	(excluded in shared counts).
	Output: aECPVlsAVlsB is an array with the ECP values for all 
	combinations of array A and B.
	NOTE: The subarrays in arrays A and B must have the same dimensions 
	(i.e. all the miRNA arrays must have the same size.).
	NOTE: Null values shall be numpy.nan.
	NOTE: aECPVlsAVlsB has arrays in A as rows and in B as columns.
	NOTE: if aALenA and aALenB ECP density is going to be calculated.
	NOTE: if aMskRef miRNA is going to be masked.
	"""
	def mltECPclc(qInJobs,qOutRslts,mthdECPclc,aMrnVlsDtA,aMrnVlsDtB, \
		fldrOutECPPrws,aALenA,aALenB,aMskRef,pntrCnts):
		"""
		Input: qInJobs is a queue with pairs of intervals. qOutRslts is
		the queue to store position in arrayA, position in arrayB, and 
		ECP value. mthdECPclc is the method to calculate the ECP value. 
		aMrnVlsDtA is an array A of arrays with values for miRNAs. 
		aMrnVlsDtB is an array B of arrays with values for miRNAs. 
		fldrOutECPPrws is a folder to store partial ECP results. aALenA 
		is an array of object lengths in the same order that aAVlsA. 
		aALenB is an array of object lengths in the same order that 
		aAVlsB. aMskRef is a mask array for the miRNAs (i.e. arrays 
		within array A and B). If pntrCnts is True  aAVlsA and aAVlsB 
		are counts so 0 values shall be considered (excluded in shared 
		counts).
		Output: qOutRslts is the queue to store position in arrayA, 
		position in arrayB, and ECP values.
		"""
		for intrvlA,intrvB in iter(qInJobs.get,'STOP'):
			lECPVlsAVlsB = mthdECPclc(aMrnVlsDtA,aMrnVlsDtB, \
			fldrOutECPPrws,intrvlA,intrvB,pntrCnts,aMskRef,aALenA, \
			aALenB)
			qOutRslts.put(lECPVlsAVlsB)
	#--------------------------
	#~ Check if there is mask for miRNAs
	if dir(aMskRef)[0]=='T':
		assert len(aMskRef) == len(aALenB[0]) == len(aALenA[0])
		if dir(aALenB)[0]=='T':
			assert dir(aALenB)[1]=='T'
			mthdECPclc = rtrnECPDnstyMskd
		else:
			assert not aALenA and not aALenB
			mthdECPclc = rtrnECPMskd
	else:
		if dir(aALenB)[0]=='T':
			assert dir(aALenB)[1]=='T'
			mthdECPclc = rtrnECPDnsty
		else:
			assert not aALenA and not aALenB
			mthdECPclc = rtrnECP
	#--------------------------
	#~ Create list of intervals for multithreading
	lenaMrnVlsDtA = len(aMrnVlsDtA)
	lenaMrnVlsDtB = len(aMrnVlsDtB)
	intrvlsMrnVlsA = []
	for strt in xrange(0,lenaMrnVlsDtA,intrvlSz):
		cEnd = strt+intrvlSz
		if cEnd<lenaMrnVlsDtA:
			end = cEnd
		else:
			end = lenaMrnVlsDtA
		intrvlsMrnVlsA.append([strt,end])
	intrvlsMrnVlsB = []
	for strt in xrange(0,lenaMrnVlsDtB,intrvlSz):
		cEnd = strt+intrvlSz
		if cEnd<lenaMrnVlsDtB:
			end = cEnd
		else:
			end = lenaMrnVlsDtB
		intrvlsMrnVlsB.append([strt,end])
	#--------------------------
	#~ Run in parallel.	
	aECPVlsAVlsB = zeros((lenaMrnVlsDtA,lenaMrnVlsDtB),dtype=float32)
	aECPVlsAVlsB.fill(nan)#fill all ECP with nan to start
	qInJobs = Queue()
	qOutRslts = Queue()
	cntVlABPrs=0	
	for intrvlA,intrvB in product(intrvlsMrnVlsA,intrvlsMrnVlsB):	
		qInJobs.put((intrvlA,intrvB))
		cntVlABPrs += 1
	for t in xrange(nThrds):
		Process(target = mltECPclc,args=(qInJobs,qOutRslts,mthdECPclc, \
		aMrnVlsDtA,aMrnVlsDtB,fldrOutECPPrws,aALenA,aALenB, \
		aMskRef,pntrCnts)).start()
	lECPVlsAVlsBGlbl = []#store global results
	for cnt in range(cntVlABPrs):
		if cnt%50==0:
			print 'Running calculations on pair %s out of %s'%(cnt, \
			cntVlABPrs)
		lECPVlsAVlsB = qOutRslts.get()
		lECPVlsAVlsBGlbl.extend(lECPVlsAVlsB)
	for t in xrange(nThrds):
		qInJobs.put('STOP')
	#--------------------------
	#~ create array: aMrnVlsDtA in rows, aMrnVlsDtB in columns.
	for vlsAPos,vlsBPos,ECP in lECPVlsAVlsBGlbl:
		aECPVlsAVlsB[vlsAPos,vlsBPos] = ECP
	if sqlFl:
		mkSqlFlECP(lECPVlsAVlsBGlbl,sqlFl,aANmsA,aANmsB)
	return aECPVlsAVlsB


########################################################
#~ Make a sqlite3 database for ECP values between genes/lncRNAs of 
# interest.
########################################################
def mkSqlFlECP(lECPVlsAVlsBGlbl,sqlFl,aANmsA,aANmsB):
	"""
	Input: lECPVlsAVlsBGlbl is a list of tuples (vrblAPos,vrblBPos,ECP). 
	vrblAPos is the position of the first variables, vrblBPos is the 
	position of the second variable, ECP is the ECP value between 
	vrblAPos and vrblBPos. A sqlite3 database will be created for the 
	input list. aANmsA is the array of variable names in the same 
	position as the numbers in vrblAPos. aANmsB is the array of variable 
	names in the same order as vrblBPos.
	Output: A sqlite3 database will be created for the input list in the 
	file sqlFl.
	"""
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	c.execute \
	('''CREATE TABLE records (id TEXT, vrblANm TEXT, vrblBNm TEXT, ECP REAL)''')
	lCnt = 0
	for vrblAPos,vrblBPos,ECP in lECPVlsAVlsBGlbl:
		vrblANm,vrblBNm = aANmsA[vrblAPos],aANmsB[vrblBPos]
		lCnt+=1
		c.execute('insert into records VALUES (?,?,?,?)', (str(lCnt), \
		vrblANm,vrblBNm,float64(ECP)))
	# create indexes. Decrease complexity of querying
	c.execute("CREATE INDEX index_records on records (id);")
	conn.commit()
	conn.close()
	return 0	


########################################################
#~ Read a sqlite3 database for correlations between genes/lncRNAs of 
# interest.
########################################################
def rtrnSqlFlECP(sqlFl,srtdVrblANms,srtdVrblBNms,rtrnECPSgnd=False):
	"""
	Input: sqlFl is a sqlite3 database with the fields id, vrblANm, 
	vrblBNm, and ECP. srtdVrblANms is a sorted lists of names 
	present in the field vrblANm. srtdVrblBNms is a sorted lists of 
	names present in the field vrblBNm. Optionally, rtrnECPSgnd can have 
	values 'negative' or 'positive', in those cases only 'negative'	or 
	'positive' ECP values are going to be retrieved respectively. 
	Output: aECPVlsAVlsB is an array of size len(srtdVrblANms) x 
	len(srtdVrblBNms) with correlation values ECP. In case the value is 
	not present nan is going to be incldued in the cell.
	NOTE: If a name is not present in a database, nan values are going 
	to be returned.
	NOTE: srtdVrblANms are going to be in rows, and srtdVrblBNms in 
	columns.
	"""
	if rtrnECPSgnd:
		try:
			if rtrnECPSgnd not in {'negative','positive'}:
				raise exceptions.CelleryWarningObjct \
				('"negative" or "positive" are values, not recognized', \
				rtrnECPSgnd)
		except exceptions.CelleryWarningObjct as err:
			print err
	#--------------------------
	#~ make a dictionary of names and positions
	lenaAVlsA = len(srtdVrblANms)
	lenaAVlsB = len(srtdVrblBNms)
	dVrblANmPos = dict([(vrblANm,pos) for pos,vrblANm in \
	enumerate(srtdVrblANms)])
	dVrblBNmPos = dict([(vrblBNm,pos) for pos,vrblBNm in \
	enumerate(srtdVrblBNms)])
	#--------------------------
	#~ make a output array
	aECPVlsAVlsB = zeros((lenaAVlsA,lenaAVlsB),dtype=float32)
	aECPVlsAVlsB.fill(nan)#fill all correlations with nan to start
	#--------------------------
	#~ retrieve variable names
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	sVrblANmsInSql = set([str(vrblANm[0]) for vrblANm in \
	c.execute('SELECT vrblANm FROM records')])
	sVrblBNmsInSql = set([str(vrblBNm[0]) for vrblBNm in \
	c.execute('SELECT vrblBNm FROM records')])
	lVrblANmInSql = list(set(srtdVrblANms).intersection(sVrblANmsInSql))
	lVrblBNmInSql = list(set(srtdVrblBNms).intersection(sVrblBNmsInSql))
	try:
		lenSrtdVrblANms = len(srtdVrblANms)
		lenlVrblANmInSql = len(lVrblANmInSql)
		if lenSrtdVrblANms!=lenlVrblANmInSql:
			raise exceptions.CelleryWarningObjct \
			('Expression for %s variable A names were retrieved out of'% \
			lenlVrblANmInSql,lenSrtdVrblANms)
	except exceptions.CelleryWarningObjct as err:
		 print err
		 pass
	try:
		lenSrtdVrblBNms = len(srtdVrblBNms)
		lenlVrblBNmInSql = len(lVrblBNmInSql)
		if lenSrtdVrblBNms!=lenlVrblBNmInSql:
			raise exceptions.CelleryWarningObjct \
			('Expression for %s variable B names were retrieved out of'% \
			lenlVrblBNmInSql,lenSrtdVrblBNms)
	except exceptions.CelleryWarningObjct as err:
		 print err
		 pass
	#--------------------------
	#~ retrieve data
	if rtrnECPSgnd == 'negative':
		cmmnd = \
		'SELECT * FROM records WHERE ECP<0 AND vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	elif rtrnECPSgnd == 'positive':
		cmmnd = \
		'SELECT * FROM records WHERE ECP>0 AND vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	else:	
		cmmnd = 'SELECT * FROM records WHERE vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	for idX,vrblANm,vrblBNm,ECP in c.execute(cmmnd):
		ECP = float32(ECP)
		aECPVlsAVlsB[dVrblANmPos[vrblANm],dVrblBNmPos[vrblBNm]] = ECP
	conn.close()
	return aECPVlsAVlsB

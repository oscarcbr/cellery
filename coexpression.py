#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  coexpression.py part of cellery (ceRNAs linking inference)
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
Methods to calculated and manipulate correlations.
"""


########################################################
#~ Import libraries.
########################################################
#----------------------------
#Import matplotlib libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#----------------------------
#Import other libraries
from cellery import exceptions
from cspearman import sprmnCrltnintrvls
from itertools import product
from multiprocessing import Queue,Process
from numpy import array,float32,float64,isnan,log2,ma,nan,zeros
from scipy.stats import normaltest,probplot

import gzip
import os
import sqlite3


########################################################
#~ Compute Spearman's correlations for all combinations of two arrays of 
# arrays with values.
########################################################
def cmpSprmCrrltn(aAVlsA,aAVlsB,aANmsA,aANmsB,tmpFldr,nThrds=10, \
	intrvlSz=700,cntRplcts=3,sqlFl=False,lMaxInQueue=4000):
	"""
	Input: aAVlsA is the array A of arrays with values. aAVlsB is the 
	array B of arrays with values. aANmsA is the array of variable names
	in the same order as aAVlsA. aANmsB is the array of variable names
	in the same order as aAVlsB. tmpFldr is a folder to store partial 
	correlation results. Optionally, nThrds is the number of threads to 
	run in parallel. Optionally, intrvlSz is the size of the interval to 
	run in multithread. Optionally, cntRplcts is the minimal number of 
	time points/replicates (or individuals) to calculate the 
	correlation. Optionally, sqlFl is a sql database to save the 
	correlations.
	Output: aCrrltnVlsAVlsB is an array with the Spearman's correlation
	in the values for all combinations of array A and B.
	NOTE: The subarrays in arrays A and B must have the same dimensions 
	(i.e. time points/individuals).
	NOTE: Null values shall be numpy.nan.
	NOTE: aCrrltnVlsAVlsB has arrays in A as rows and in B as columns.
	"""
	def mltSprmn(qInJobs,qOutRslts,aAVlsA,aAVlsB,tmpFldr,cntRplcts):
		"""
		Input: qInJobs is a queue with pairs of intervals. qOutRslts is
		the queue to store position in arrayA, position in arrayB, and 
		Spearman's correlation. aAVlsA is the array A of arrays with 
		values. aAVlsB is the array B of arrays with values. tmpFldr is 
		a folder to store partial correlation results. cntRplcts is the 
		minimal number of time points/replicates to calculate the 
		correlation.
		Output: qOutRslts is the queue to store position in arrayA, 
		position in arrayB, and Spearman's correlation.
		"""
		for intrvlA,intrvB in iter(qInJobs.get,'STOP'):
			lSprmnVlsAVlsB = sprmnCrltnintrvls(aAVlsA,aAVlsB, \
			intrvlA,intrvB,tmpFldr,cntRplcts)
			qOutRslts.put(lSprmnVlsAVlsB)
	#--------------------------
	#~ Create list of intervals for multithreading
	lenaAVlsA = len(aAVlsA)
	lenaAVlsB = len(aAVlsB)
	intrvlsVlsA = []
	for strt in xrange(0,lenaAVlsA,intrvlSz):
		cEnd = strt+intrvlSz
		if cEnd<lenaAVlsA:
			end = cEnd
		else:
			end = lenaAVlsA
		intrvlsVlsA.append([strt,end])
	intrvlsVlsB = []
	for strt in xrange(0,lenaAVlsB,intrvlSz):
		cEnd = strt+intrvlSz
		if cEnd<lenaAVlsB:
			end = cEnd
		else:
			end = lenaAVlsB
		intrvlsVlsB.append([strt,end])	
	#--------------------------
	#~ Set parameteres to run in parallel
	lQueueIn = []#holder for queues
	aCrrltnVlsAVlsB = zeros((lenaAVlsA,lenaAVlsB),dtype=float32)
	aCrrltnVlsAVlsB.fill(nan)#fill all correlations with nan to start
	qInJobs = Queue()
	cntVlABPrs=0
	#--------------------------
	#~ Limit the size of input queues to lMaxInQueue
	for intrvlA,intrvB in product(intrvlsVlsA,intrvlsVlsB):	
		qInJobs.put((intrvlA,intrvB))
		cntVlABPrs += 1
		if cntVlABPrs>=lMaxInQueue:
			lQueueIn.append(qInJobs)
			qInJobs = Queue()
			cntVlABPrs=0
	if cntVlABPrs<lMaxInQueue:
		lQueueIn.append(qInJobs)
	#--------------------------
	#~ Run in parallel
	lenlQueueIn = len(lQueueIn)+1
	while lQueueIn:
		lenlQueueIn-=1
		qOutRslts = Queue()
		qInJobs = lQueueIn.pop()
		cntVlABPrs = int(qInJobs.qsize())
		for t in xrange(nThrds):
			Process(target = mltSprmn,args=(qInJobs,qOutRslts,aAVlsA, \
			aAVlsB,tmpFldr,cntRplcts)).start()
		for cnt in xrange(cntVlABPrs):
			if cnt%50==0:
				print \
				'Running calculations on pair %s out of %s, of queue %s' \
				%(cnt,cntVlABPrs,lenlQueueIn)
			lSprmnVlsAVlsB = qOutRslts.get()
			#--------------------------
			#~ create array: aAVlsA in rows, aAVlsB in columns.
			for vlsAPos,vlsBPos,sprmnCrltn in lSprmnVlsAVlsB:
				aCrrltnVlsAVlsB[vlsAPos,vlsBPos] = sprmnCrltn
		for t in xrange(nThrds):
			qInJobs.put('STOP')
	#--------------------------
	#~ create sql database
	if sqlFl:
		mkSqlFlCrrltnFrmArry(aCrrltnVlsAVlsB,sqlFl,aANmsA,aANmsB)
	return aCrrltnVlsAVlsB


########################################################
#~ Read expression values from a sqlite3 table for genes/lncRNAs or 
# regions of interest and different times/individuals.
########################################################
def cmpCrltnSqlFl(inSqlFl,sGnIntrst,outSqlFl=False):
	"""
	Input: inSqlFl is a sqlite3 database with the attributes of interest.
	sGnIntrst is a set with names of gene identifiers/names of interest.
	Output: dNmTrpleMrgdChrLIntrvlsStrndCmmnm is a dictionary of 
	gene/transcript identifier/name in column 0 -or posNm- as keys and 
	as values tuples of (chromosome, list of merged intervals, strand, 
	and common names).
	NOTE: starts are 0-based and ends 1-based. 
	NOTE: id is the field for gene identifiers/names.
	"""
	dNmTrpleMrgdChrLIntrvlsStrndCmmnm = {}
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	for info in c.execute('SELECT * FROM records'):
		gn,chr,strts,ends,strnd,cmmNm = [str(clmn) for clmn in info]
		strts = [int(s) for s in strts.split(';')]
		ends = [int(e) for e in ends.split(';')]
		aMrdgIntrvls = zip(strts,ends)
		dNmTrpleMrgdChrLIntrvlsStrndCmmnm[gn] = (chr,aMrdgIntrvls, \
		strnd,cmmNm)
	if sGnIntrst:
		sGnDB = set(dNmTrpleMrgdChrLIntrvlsStrndCmmnm.keys())
		sGnsDiff = sGnIntrst.difference(sGnDB)
		sGnsRmv = sGnDB.difference(sGnIntrst)
		try:
			if sGnsDiff:
				raise exceptions.CelleryWarningObjct \
				('%s identifiers not present in'%len(sGnsDiff),sqlFl)
		except exceptions.CelleryWarningObjct as err:
			 print err
			 pass
		while sGnsRmv:
			gn = sGnsRmv.pop()
			rmvd = dNmTrpleMrgdChrLIntrvlsStrndCmmnm.pop(gn)
	return dNmTrpleMrgdChrLIntrvlsStrndCmmnm


########################################################
#~ Calculate correlations for genes, taking as replicates the expression
# of invdividuals. Takes as input GTeX results and calculates 
# correlations for each tissue.
########################################################
def grpGTxSmplByTiss(fExprssn,fGTeXSmplAttrbt,fGTeXSbjctAttrbt, \
	sGnIntrst,GTeXCrrltnFldr,frzn=True,frznClmn=17,tissTypClmn=6, \
	sqlDBextsn='.db',sTissIntrst=False,nThrds=10):
	"""
	Input: fExpressn is the full path to the file with the expression 
	from GTeX database. fGTeXSmplAttrbt is the file with the sample 
	attribute as downloaded from GTeX. fGTeXSbjctAttrbt is the file with 
	the subjects' attributes as downloaded from GTeX. sGnIntrst is a set 
	of gene or transcript of interest. GTeXCrrltnFldr is the output 
	folder where the results are going to be written. Optionally, if 
	frzn is True, only the samples used in the frozen version of GTeX 
	are going to be included. frznClmn is the column with the value for 
	the frozen status. tissTypClmn is the column with the tissue type. 
	sqlDBextsn is the extension to save the databases with the results 
	from each tissue. sTissIntrst is a set of tissues of interest, if 
	not false only results for those tissues are going to be retrieved.
	nThrds is the number of threads to run in parallel the calculations 
	of the correlations.
	Output: Databases with the correlations for each pair of genes for 
	each tissue of interest. Additionally, a dictionary with tissues as
	keys and array of gene correlations (sorted by names) is retrieved 
	(i.e. dTissACrrltnGnsTiss).
	NOTE: If frzn is True, only the samples used in the frozen version
	of GTeX are going to be included. These samples have been identified 
	as those best suited for use in analysis with a specific focus on 
	eQTL analysis.
	NOTE: The minimum number fo samples individual required is 3 by 
	default. This parameter can be changed in cntRplcts method after 
	deleting the databases.
	NOTE: The expression of each sample is log2-normalized.
	"""
	#--------------------------
	#~ Make a dictionary and array of individuals
	lIndvls = []
	hdr = True
	for el in open(fGTeXSbjctAttrbt,'r'):
		if hdr:
			hdr = False
		elif el.strip():
			sbjctId = el.split('\t')[0]
			lIndvls.append(sbjctId)		
	dSbjctPos = dict([(indvl,pos) for pos,indvl in enumerate(lIndvls)])
	#--------------------------
	#~ Make a dictionary of tissues and individuals
	dSmplIndvl,dSmplTiss = {},{}
	hdr = True
	for el in open(fGTeXSmplAttrbt,'r'):
		if hdr:
			hdr = False
		elif el.strip():
			inf=el.split('\t')
			#decide if include the sample
			if frzn:
				if inf[17]=='USE ME':
					appndSmpl = True
				else:
					appndSmpl = False
			else:
				appndSmpl = True
			if appndSmpl:
				smpl = inf[0]
				indvl = '-'.join(smpl.split('-')[:2])
				tiss = inf[tissTypClmn]
				dSmplIndvl[smpl] = indvl
				if sTissIntrst:
					if tiss in sTissIntrst:
						dSmplTiss[smpl] = tiss
				else:
					dSmplTiss[smpl] = tiss
	srtdTiss = sorted(set(dSmplTiss.values()))
	if sTissIntrst:
		assert set(srtdTiss)==sTissIntrst
	#--------------------------
	#~ Sort individuals, tissues and genes in each sample
	srtGnNms = sorted(sGnIntrst)
	dGnNmPos = dict([(gn,pos) for pos,gn in enumerate(srtGnNms)])
	nGns = len(sGnIntrst)
	nIndvls = len(lIndvls)
	nTiss = len(srtdTiss)
	dTissGnsIndvls =  dict([(tiss,zeros((nGns,nIndvls),dtype=float32)) \
	for tiss in xrange(nTiss)])
	dTissIndvlsCnts = dict([(tiss,zeros(nIndvls,dtype=float32)) for tiss \
	in xrange(nTiss)])
	#--------------------------
	#~ Fill the array and calculate mean expression
	hdr = True
	for el in gzip.open(fExprssn,'rb'):
		if el.strip():
			inf=el.splitlines()[0].split()
			if hdr:
				hdr = False
				infSmpls = inf[4:]
				lenData = len(infSmpls)
				aTissInArry,aIndvlInArry = [],[]
				for smpl in infSmpls:
					posInIndvl = lIndvls.index(dSmplIndvl[smpl])
					if sTissIntrst:
						if dSmplTiss.has_key(smpl):
							posInTiss = srtdTiss.index(dSmplTiss[smpl])
							dTissIndvlsCnts[posInTiss][posInIndvl]+=1
						else:
							posInTiss = nan
					else:
						posInTiss = srtdTiss.index(dSmplTiss[smpl])
						dTissIndvlsCnts[posInTiss][posInIndvl]+=1
					aTissInArry.append(posInTiss)
					aIndvlInArry.append(posInIndvl)
				aTissInArry = array(aTissInArry)
				aIndvlInArry = array(aIndvlInArry)
			else:
				gnNm = inf[0].split('.')[0]				
				if gnNm in sGnIntrst:
					gnPos = dGnNmPos[gnNm]
					aExprsn = array(inf[4:],dtype=float32)
					for pos in xrange(lenData):
						posInTiss = aTissInArry[pos]
						if not isnan(posInTiss):
							exprsn = aExprsn[pos]
							dTissGnsIndvls[posInTiss][gnPos] \
							[aIndvlInArry[pos]] += exprsn
	#--------------------------
	#~ Calculate correlations (of log2 expressions) and make database
	dTissACrrltnGnsTiss = {}#output dictionary
	for pos,tiss in enumerate(srtdTiss):
		print '--------------------------------------------------------'
		print 'Running calculations for tissue %s, number %s out of %s'% \
		(tiss,pos,nTiss)
		appndTiss = False
		if sTissIntrst:
			if tiss in sTissIntrst:
				appndTiss = True
		else:
			appndTiss = True
		if appndTiss:
			outGTeXDBFl = os.path.join(GTeXCrrltnFldr,''.join([tiss. \
			replace(' ','_'),sqlDBextsn]))
			outGTeXDBNrmPlt = os.path.join(GTeXCrrltnFldr,''.join([tiss. \
			replace(' ','_'),'.png']))
			if not os.path.exists(outGTeXDBFl):
				aGnsIndvls = dTissGnsIndvls.pop(pos)
				aIndvlsCnts = dTissIndvlsCnts.pop(pos)
				print aGnsIndvls
				print aIndvlsCnts
				aGnsIndvls/=aIndvlsCnts#mean
				aGnsIndvls = log2(aGnsIndvls)
				print aGnsIndvls
				aGnsIndvls = ma.masked_invalid(aGnsIndvls)
				aGnsIndvls.fill_value = nan
				if not aGnsIndvls.mask.all():#at least there is one valid
					#--------------------------
					#~ Make normality calculations
					aGnsIndvlsCmprsd = aGnsIndvls.compressed()
					k2,pval = normaltest(aGnsIndvlsCmprsd)
					print '--------------------------------------------------------'
					print 'Normality test for tissue: %s. k2: %s and p-val: %s...'% \
					(tiss,k2,pval)
					print '--------------------------------------------------------'
					fig = plt.figure()
					ax = fig.add_subplot(111)
					probplot(aGnsIndvlsCmprsd, dist="norm",	plot=plt)
					plt.gcf().subplots_adjust(bottom=0.15)
					plt.savefig(outGTeXDBNrmPlt,bbox_inches='tight')
					plt.close()	
					del(aGnsIndvlsCmprsd)#to free memory
					#--------------------------
					aGnNms = array(srtGnNms)
					aGnsIndvls = aGnsIndvls.filled()
					aCrrltnGnsTiss = cmpSprmCrrltn(aGnsIndvls,aGnsIndvls, \
					aGnNms,aGnNms,GTeXCrrltnFldr,sqlFl=outGTeXDBFl, \
					nThrds=nThrds)
					os.system('rm %s\*.ecp'%GTeXCrrltnFldr)
				else:
					aCrrltnGnsTiss = None
			else:
				aGnNms = array(srtGnNms)
				aCrrltnGnsTiss = rtrnSqlFlCrrltn(outGTeXDBFl,aGnNms, \
				aGnNms)
			dTissACrrltnGnsTiss[tiss] = aCrrltnGnsTiss
		print '--------------------------------------------------------'
	return dTissACrrltnGnsTiss


########################################################
#~ Calculate correlations for genes, taking a replicates different 
# data points. Input is a ".tsv" file.
########################################################
def grpTsvSmplByPoint(fExprssnTsv,sGnIntrst,tsvSmplCrrltnFldr,dbNm, \
	sqlDBextsn='.db',srtPntsIntrst=False,nThrds=10):
	"""
	Input: fExpressn is the full path to the file with the normalized 
	expression in "tsv" format. sGnIntrst is a set of gene or transcript 
	of interest. tsvSmplCrrltnFldr is the full path to the output folder. 
	dbNm is the prefix for the name of the database. sqlDBextsn is the 
	extension to save the databases with the results from each tissue. 
	srtPntsIntrst is the name of the sorted data points of interest. 
	nThrds is the number of threads to run in parallel the calculations 
	of the correlations.
	Output: An array with the correlations for each pair of genes for 
	the input data. Additionally, a database is created in the output
	folder.
	NOTE: The input file must have a header.
	NOTE: If srtPntsIntrst of the columns in the input file that matches
	the expressions in srtPntsIntrst are included.
	NOTE: The minimum number fo samples individual required is 3 by 
	default. This parameter can be changed in cntRplcts method after 
	deleting the databases.
	NOTE: The expression of each sample is log2-normalized.
	"""
	#--------------------------
	#~ Sort genes of interest
	srtGnNms = sorted(sGnIntrst)
	dGnNmPos = dict([(gn,pos) for pos,gn in enumerate(srtGnNms)])
	nGns = len(sGnIntrst)
	#--------------------------
	#~ Fill the array and calculate mean expression
	hdr = True
	for el in open(fExprssn,'r'):
		if el.strip():
			inf=el.splitlines()[0].split('\t')
			if hdr:
				hdr = False
				infSmpls = inf[1:]
				lenData = len(infSmpls)
				if sPntsIntrst:
					sPntsIntrst = set(srtPntsIntrst)
					aPnts = array([srtPntsIntrst.index(smpl) for smpl \
					in infSmpls if smpl in sPntsIntrst])
				else:
					aPnts = array([pos for pos,smpl in \
					enumerate(infSmpls)])
				nPnts = len(aPnts)
				aGnsPnts = zeros((nGns,nPnts),dtype=float32)				
			else:
				gnNm = inf[0].split('.')[0]				
				if gnNm in sGnIntrst:
					gnPos = dGnNmPos[gnNm]
					aExprsn = array(inf[1:],dtype=float32)	
					for pos in xrange(lenData):
						exprsn = aExprsn[pos]
						aGnsPnts[gnPos][aPnts[pos]] = exprsn
	#--------------------------
	#~ Calculate correlations and make database
	outGTeXDBFl = os.path.join(GTeXCrrltnFldr,''.join([dbNm,sqlDBextsn]))
	outGTeXDBNrmPlt = os.path.join(GTeXCrrltnFldr,''.join([dbNm,'.png']))
	if not os.path.exists(outGTeXDBFl):
		aGnsPnts = ma.masked_invalid(aGnsPnts)
		aGnsPnts.fill_value = nan
		if not aGnsPnts.mask.all():#at least there is one valid
			#--------------------------
			#~ Make normality calculations
			aGnsPntsCmprsd = aGnsPnts.compressed()
			k2,pval = normaltest(aGnsPntsCmprsd)
			print '--------------------------------------------------------'
			print 'Normality test for tissue: %s. k2: %s and p-val: %s...'% \
			(tiss,k2,pval)
			print '--------------------------------------------------------'
			fig = plt.figure()
			ax = fig.add_subplot(111)
			probplot(aGnsPntsCmprsd, dist="norm",	plot=plt)
			plt.gcf().subplots_adjust(bottom=0.15)
			plt.savefig(outGTeXDBNrmPlt,bbox_inches='tight')
			plt.close()	
			del(aGnsPntsCmprsd)#to free memory
			#--------------------------
			aGnNms = array(srtGnNms)
			aGnsPnts = aGnsPnts.filled()
			aCrrltnGnsTiss = cmpSprmCrrltn(aGnsPnts,aGnsPnts, \
			aGnNms,aGnNms,GTeXCrrltnFldr,sqlFl=outGTeXDBFl, \
			nThrds=nThrds)
			os.system('rm %s\*.ecp'%GTeXCrrltnFldr)
		else:
			aCrrltnGnsTiss = None
	else:
		aGnNms = array(srtGnNms)
		aCrrltnGnsTiss = rtrnSqlFlCrrltn(outGTeXDBFl,aGnNms,aGnNms)
	return aCrrltnGnsTiss


########################################################
#~ Read correlations from a TSV file
########################################################
def mkCrrltnFrmTSVfl(fCrrltnTsv):
	"""
	Input: fCrrltnTsv is a file in "tsv" format with header have name 
	of genes. In the same way, column 0 shall have gene names.
	Output: aCrrltnVlsAVlsB is an array with the Spearman's correlation
	in the values for all combinations of array A and B. aANmsB is the 
	array of variable names in the same order as the header. aANmsA is 
	the array of variable names in the same order as column 0. 
	"""	
	#--------------------------
	#~ Fill the array and calculate mean expression
	hdr = True
	aANmsA = []
	aDtNmANmB = []
	for el in open(fExprssn,'r'):
		if el.strip():
			inf=el.splitlines()[0].split('\t')
			if hdr:
				hdr = False
				infSmpls = inf[1:]
				lenData = len(infSmpls)
				aANmsB = array([pos for pos,smpl in \
				enumerate(infSmpls)])
				lenDtNmsB = len(aANmsB)
			else:
				appndInf = False
				gnNm = inf[0].split('.')[0]	
				aANmsA.append(gnNm)
				aExprsn = array(inf[1:],dtype=float32)	
				for pos in xrange(lenDtNmsA):
					exprsn = aExprsn[pos]
					aDtNmANmB.append(exprsn)
	#--------------------------
	#~ Make output array
	lenDtNmsA = len(aANmsA)
	aCrrltnVlsAVlsB = zeros((lenDtNmsA,lenDtNmsB),dtype=float32)
	for posDtNmB in xrange(lenDtNmsB):
		for posDtNmA in xrange(lenDtNmsA):
			exprssVl = aDtNmANmB.pop(0)
			aCrrltnVlsAVlsB[posDtNmA][posDtNmB][posInCrrltnList] = \
			exprssVl
	return aCrrltnVlsAVlsB,aANmsA,aANmsB
		

########################################################
#~ Make a sqlite3 database for correlations between genes/lncRNAs of 
# interest.
########################################################
def mkSqlFlCrrltn(lSprmnVlsAVlsBGlbl,sqlFl,aANmsA,aANmsB):
	"""
	Input: lSprmnVlsAVlsBGlbl is a list of tuples (vrblAPos, vrblBPos,
	sprmnCrltn). vrblAPos is the position of the first variables, 
	vrblBPos is the position of the second variable, sprmnCrltn is the 
	correlation between vrblAPos and vrblBPos. A sqlite3 database will 
	be created for the input list. aANmsA is the array of variable names 
	in the same position as the numbers in vrblAPos. aANmsB is the array 
	of variable names in the same order as vrblBPos. 
	Output: A sqlite3 database will be created for the input list in the 
	file sqlFl.
	"""
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	c.execute \
	('''CREATE TABLE records (id TEXT, vrblANm TEXT, vrblBNm TEXT, sprmnCrltn REAL)''')
	records = []
	lCnt = 0
	for vrblAPos,vrblBPos,sprmnCrltn in lSprmnVlsAVlsBGlbl:
		vrblANm,vrblBNm = aANmsA[vrblAPos],aANmsB[vrblBPos]
		lCnt+=1
		records.append((str(lCnt),vrblANm,vrblBNm,float64(sprmnCrltn)))
	records = tuple(records)
	c.executemany('insert into records VALUES (?,?,?,?)', records)
	# create indexes. Decrease complexity of querying
	c.execute("CREATE INDEX index_records on records (id);")
	conn.commit()
	conn.close()
	return 0	
		

########################################################
#~ Make a sqlite3 database for correlations between genes/lncRNAs of 
# interest.
########################################################
def mkSqlFlCrrltnFrmArry(aCrrltnVlsAVlsB,sqlFl,aANmsA,aANmsB):
	"""
	Input: aCrrltnVlsAVlsB is an of correlations with size len(aANmsA)X
	len(aANmsB). A sqlite3 database will be created for the input list
	in the file sqlFl. aANmsA is the array of variable names in the same 
	position as the numbers in vrblAPos. aANmsB is the array of variable 
	names in the same order as vrblBPos. 
	Output: A sqlite3 database will be created for the input list in the 
	file sqlFl.
	"""
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	c.execute \
	('''CREATE TABLE records (id TEXT, vrblANm TEXT, vrblBNm TEXT, sprmnCrltn REAL)''')
	lCnt = 0
	for vrblAPos in xrange(len(aANmsA)):
		for vrblBPos in xrange(len(aANmsB)):
			vrblANm,vrblBNm = aANmsA[vrblAPos],aANmsB[vrblBPos]
			sprmnCrltn = aCrrltnVlsAVlsB[vrblAPos][vrblBPos]
			if not isnan(sprmnCrltn):
				lCnt+=1
				c.execute('insert into records VALUES (?,?,?,?)', (str(lCnt), \
				vrblANm,vrblBNm,float64(sprmnCrltn)))
	# create indexes. Decrease complexity of querying
	c.execute("CREATE INDEX index_records on records (id);")
	conn.commit()
	conn.close()
	return 0	


########################################################
#~ Read a sqlite3 database for correlations between genes/lncRNAs of 
# interest.
########################################################
def rtrnSqlFlCrrltn(sqlFl,srtdVrblANms,srtdVrblBNms, \
	rtrnSgndCrrltn=False):
	"""
	Input: sqlFl is a sqlite3 database with the fields id, vrblANm, 
	vrblBNm, and sprmnCrltn. srtdVrblANms is a sorted lists of names 
	present in the field vrblANm. srtdVrblBNms is a sorted lists of 
	names present in the field vrblBNm. Optionally, rtrnSgndCrrltn can
	have values 'negative' or 'positive', in those cases only 'negative'
	or 'positive' correlations are going to be retrieved respectively.
	Output: aCrrltnVlsAVlsB is an array of size len(srtdVrblANms) x 
	len(srtdVrblBNms) with correlation values sprmnCrltn. In case the 
	value is not present nan is going to be incldued in the cell.
	NOTE: If names is not present in a database, a nan value is going to
	be returned.
	NOTE: srtdVrblANms are going to be in rows, and srtdVrblBNms in 
	columns.
	"""
	if rtrnSgndCrrltn:
		try:
			if rtrnSgndCrrltn not in {'negative','positive'}:
				raise exceptions.CelleryWarningObjct \
				('"negative" or "positive" are values, not recognized', \
				rtrnSgndCrrltn)
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
	aCrrltnVlsAVlsB = zeros((lenaAVlsA,lenaAVlsB),dtype=float32)
	aCrrltnVlsAVlsB.fill(nan)#fill all correlations with nan to start
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
	if rtrnSgndCrrltn == 'negative':
		cmmnd = \
		'SELECT * FROM records WHERE sprmnCrltn<0 AND vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	elif rtrnSgndCrrltn == 'positive':
		cmmnd = \
		'SELECT * FROM records WHERE sprmnCrltn>0 AND vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	else:	
		cmmnd = 'SELECT * FROM records WHERE vrblANm in (%s) AND vrblBNm in (%s)' \
		%(','.join(['"%s"'%v for v in lVrblANmInSql]),','.join(['"%s"'%v \
		for v in lVrblBNmInSql]))
	for idX,vrblANm,vrblBNm,sprmnCrltn in c.execute(cmmnd):
		sprmnCrltn = float32(sprmnCrltn)
		aCrrltnVlsAVlsB[dVrblANmPos[vrblANm],dVrblBNmPos[vrblBNm]] = \
		sprmnCrltn
	conn.close()
	return aCrrltnVlsAVlsB



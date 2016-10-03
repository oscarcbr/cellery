#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parser.py part of cellery (ceRNAs linking inference)
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
Libraries with parasers for different MRE predicting algorithms for cellery
"""


########################################################
#~ Import external libraries
########################################################
from cellery import exceptions
from numpy import nan,float32,zeros

import os


########################################################
#~ Parse mirMap results
########################################################
def prsMirMapRslts(mirMapRsltsFldr,srtdMirnaFms,extnsnMirMap,lDtObjct):
	"""
	Input: mirMapRsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted 
	list of miRNAs of interest. extnsnMirMap is the extension of the 
	files in the results folder. lDtObjct is a list of data objects, 
	must be the same as that of the files in mirMapRsltsFldr excepting 
	for the extension. The objects in this list must include the pointer 
	"pos", "name" and "aMirnaNms" (in a congruent way with the input 
	list). nThrds is the number of threads to run in multithread. 
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", "aMirnaMirMapScrs" (in the same order as 
	"aMirnaNms") for the genes/lncrnas of interest in mirMapRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsMirMapFile(mirMapRsltsFl,srtdMirnaFms):
		"""
		Input: "mirMapRsltsFl" is an output file from mirMap program. 
		"srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaMirMapCnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". aMirnaMirMapScrs are the scores for MRES of 
		miRNAs in "srtdMirnaFms". They are in the same order as in 
		"srtdMirnaFms".
		"""
		dMirMapScrs = {}
		dMirMapCnts = {}
		strtng = True
		cnt = 0
		for eline in open(mirMapRsltsFl,'r'):
			if eline.find('>')==0:
				if strtng:
					strtng = False
				else:
					dMirMapScrs[mirna] = scr
					dMirMapCnts[mirna] = cnt
					cnt = 0
				mirna = eline[1:].splitlines()[0].replace('_','-')
				scr = nan
			elif eline.find('mirMap_Score')==0:
				scr = float32(eline.split('mirMap_Score:\t')[1]. \
				splitlines()[0].strip())
			elif eline.find('UTR position')>-1:
				cnt+=1
		#add last case
		if not strtng:
			dMirMapScrs[mirna] = scr	
			dMirMapCnts[mirna] = cnt	
		#create aMirna object
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaMirMapScrs = zeros(cntMirnaNms,dtype = float32)
		aMirnaMirMapCnts = zeros(cntMirnaNms,dtype = float32)
		aMirnaMirMapScrs.fill(nan)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirMapScrs.has_key(mirna):
				aMirnaMirMapScrs[mirPos] = dMirMapScrs[mirna]
				aMirnaMirMapCnts[mirPos] = dMirMapCnts[mirna]
		return aMirnaMirMapScrs,aMirnaMirMapCnts
	#----------------------------
	#Parse mirMap file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		mirMapRsltsFl = os.path.join(mirMapRsltsFldr, \
		''.join([dt.name,extnsnMirMap]))
		if os.path.exists(mirMapRsltsFl):
			aMirnaMirMapScrs,aMirnaMirMapCnts = \
			prsMirMapFile(mirMapRsltsFl,srtdMirnaFms)
		else:
			aMirnaMirMapScrs = zeros(cntMirnaNms,dtype = float32)
			aMirnaMirMapCnts = zeros(cntMirnaNms,dtype = float32)
			aMirnaMirMapScrs.fill(nan)
		dt.aMirnaMirMapScrs = aMirnaMirMapScrs	
		dt.aMirnaMirMapCnts = aMirnaMirMapCnts
	return 0


########################################################
#~ Parse miranda results
########################################################
def prsMirandaRslts(mirandaRsltsFldr,srtdMirnaFms,extnsnMiranda, \
	lDtObjct):
	"""
	Input: mirandaRsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted list of 
	miRNAs of interest. extnsnMiranda is the extension of the files in 
	the results folder. lDtObjct is a list of data objects, must be 
	the same as that of the files in mirandaRsltsFldr excepting for the 
	extension. The objects in this list must include the pointer "pos", 
	"name" and "aMirnaNms" (in a congruent way with the input list).
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", "aMirnaMirndCnts", "aMirnaMirndScrs", 
	and "aMirnaMirndEngy" (in the same order as "aMirnaNms") for the 
	genes/lncrnas of interest in mirandaRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsMirandaFile(mirandaRsltsFl,srtdMirnaFms):
		"""
		Input: "mirandaRsltsFl" is an output file from miranda program. 
		"srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaMirndCnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". aMirnaMirndEngy are energies of non-overlapping 
		MREs for miRNAs in "srtdMirnaFms". aMirnaMirndScrs are alignment 
		scores of MREs for miRNAs in "srtdMirnaFms". They are in the 
		same order as in "srtdMirnaFms".
		"""
		dMirnasMirandaScrs = {}
		dMirnasMirandaEnrgs = {}
		dMirnasMirandaCnts = {}
		#there is synchronization problem that I haven't been able to 
		#resolve excepting by this mean... no documentation on python 
		#about it
		trls = 30
		while trls > 0:
			try:
				openfile = open(mirandaRsltsFl,'r')
				break
			except:
				trls -= 1
				print '\t...Waiting for the file to update...'
				sleep(2)
		if trls == 0:
			 raise exceptions.CelleryExceptionObjct \
			 ('\nCould not be read the file %s!'%mirandaRsltsFl)
		for eline in openfile:
			if eline.strip() and eline[0] == '>' and eline[1] == '>':
				allinf = eline.split('\t')
				mirnNm = allinf[0]
				mirna = '-'.join(mirnNm.split('-')[1:])
				maxScore = float32(allinf[4])
				minFrEnrgy = float32(allinf[5])
				cntMtchs = len(allinf[9].split())
				dMirnasMirandaScrs[mirna] = maxScore
				dMirnasMirandaEnrgs[mirna] = minFrEnrgy
				dMirnasMirandaCnts[mirna] = cntMtchs
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaMirndCnts = zeros(cntMirnaNms,dtype = float32)
		aMirnaMirndScrs = zeros(cntMirnaNms,dtype = float32)
		aMirnaMirndEngy = zeros(cntMirnaNms,dtype = float32)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirnasMirandaCnts.has_key(mirna):
				cntMtchs = dMirnasMirandaCnts.pop(mirna)
				minFrEnrgy = dMirnasMirandaEnrgs.pop(mirna)
				maxScore = dMirnasMirandaScrs.pop(mirna)
				aMirnaMirndCnts[mirPos] = cntMtchs
				aMirnaMirndEngy[mirPos] = minFrEnrgy
				aMirnaMirndScrs[mirPos] = maxScore
			else:
				aMirnaMirndEngy[mirPos] = nan
				aMirnaMirndScrs[mirPos] = nan
		return aMirnaMirndCnts,aMirnaMirndScrs,aMirnaMirndEngy
	#----------------------------
	#Parse Miranda file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		mirandaRsltsFl = os.path.join(mirandaRsltsFldr, \
		''.join([dt.name,extnsnMiranda]))
		if os.path.exists(mirandaRsltsFl):
			aMirnaMirndCnts,aMirnaMirndScrs,aMirnaMirndEngy = \
			prsMirandaFile(mirandaRsltsFl,srtdMirnaFms)
		else:
			aMirnaMirndCnts = zeros(cntMirnaNms,dtype = float32)
			aMirnaMirndScrs = zeros(cntMirnaNms,dtype = float32)
			aMirnaMirndEngy = zeros(cntMirnaNms,dtype = float32)
			aMirnaMirndScrs.fill(nan)
			aMirnaMirndEngy.fill(nan)
		dt.aMirnaMirndCnts = aMirnaMirndCnts
		dt.aMirnaMirndScrs = aMirnaMirndScrs
		dt.aMirnaMirndEngy = aMirnaMirndEngy
	return 0


########################################################
#~ Parse PITA results
########################################################
def prsPITARslts(PITARsltsFldr,srtdMirnaFms,extnsnPITA,lDtObjct):
	"""
	Input: PITARsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted list 
	of miRNAs of interest. extnsnPITA is the extension of the files in 
	the results folder. lDtObjct is a list of data objects, must be 
	the same as that of the files in PITARsltsFldr excepting for the 
	extension. The objects in this list must include the pointer "pos", 
	"name" and "aMirnaNms" (in a congruent way with the input list).
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", "aMirnaPITAScrs" and "aMirnaPITACnts" 
	(in the same order as "aMirnaNms") for the genes/lncrnas of interest 
	in PITARsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsPITAFile(PITARsltsFl,srtdMirnaFms):
		"""
		Input: "PITARsltsFl" is an output file from PITA program. 
		"srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaPITACnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". aMirnaPITAScrs are the scores for MREs in 
		"srtdMirnaFms". They are in the same order as in "srtdMirnaFms".
		"""
		dMirnasPITACnts = {}
		dMirnasPITAScrs = {}
		hdr = True
		for eline in open(PITARsltsFl,'r'):
			if hdr:
				hdr = False
			elif eline.strip():
				eline = eline.strip()
				allinf = eline.split('\t')
				mirnNm = allinf[1].split()[0]
				mirna = '-'.join(mirnNm.split('-')[1:])
				cntMtchs = int(allinf[2])
				scr = float32(allinf[3])
				dMirnasPITACnts[mirna] = cntMtchs
				dMirnasPITAScrs[mirna] = scr
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaPITACnts = zeros(cntMirnaNms,dtype = float32)
		aMirnaPITAScrs = zeros(cntMirnaNms,dtype = float32)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirnasPITACnts.has_key(mirna):
				cntMtchs = dMirnasPITACnts.pop(mirna)
				aMirnaPITACnts[mirPos] = cntMtchs
				scr = dMirnasPITAScrs.pop(mirna)
				if scr>0:
					aMirnaPITAScrs[mirPos] = scr
				else:
					aMirnaPITAScrs[mirPos] = nan
			else:
				aMirnaPITAScrs[mirPos] = nan
		return aMirnaPITACnts,aMirnaPITAScrs
	#----------------------------
	#Parse PITA file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		PITARsltsFl = os.path.join(PITARsltsFldr, \
		''.join([dt.name,extnsnPITA]))
		if os.path.exists(PITARsltsFl):
			aMirnaPITACnts,aMirnaPITAScrs =	prsPITAFile(PITARsltsFl, \
			srtdMirnaFms)
		else:
			aMirnaPITACnts = zeros(cntMirnaNms,dtype = float32)
			aMirnaPITAScrs = zeros(cntMirnaNms,dtype = float32)
			aMirnaPITAScrs.fill(nan)
		dt.aMirnaPITACnts = aMirnaPITACnts
		dt.aMirnaPITAScrs = aMirnaPITAScrs
	return 0


########################################################
#~ Parse RNAhybrid results
########################################################
def prsRNAhybrdRslts(RNAhybrdRsltsFldr,srtdMirnaFms,extnsnRNAhybrd, \
	lDtObjct):
	"""
	Input: RNAhybrdRsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted list of 
	miRNAs of interest. extnsnRNAhybrd is the extension of the files in 
	the results folder. lDtObjct is a list of data objects, shall be 
	the same as that of the files in RNAhybrdRsltsFldr excepting for the 
	extension. The objects in this list must include the pointer "pos", 
	"name" and "aMirnaNms" (in a congruent way with the input list).
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", "aMirnaRNAhCnts", and "aMirnaRNAhEngy" 
	(in the same order as "aMirnaNms") for the genes/lncrnas of interest 
	in RNAhybrdRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsRNAhybrdFile(RNAhybrdRsltsFl,srtdMirnaFms):
		"""
		Input: "RNAhybrdRsltsFl" is an output file from RNAhybrid 
		program. "srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaRNAhCnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". aMirnaRNAhEngy are energies of non-overlapping 
		MREs for miRNAs in "srtdMirnaFms". They are in the same order as 
		in "srtdMirnaFms".
		"""
		dMirnasRNAhybrdScrs = {}
		for eline in open(RNAhybrdRsltsFl,'r'):
			if eline.strip():
				allinf = eline.split(':')
				if allinf[0] != 'target too long' and len(allinf)>2:
					mirnNm = allinf[2]
					mirna = '-'.join(mirnNm.split('-')[1:])
					frEnrgy = float32(allinf[4])
					if dMirnasRNAhybrdScrs.has_key(mirna):
						dMirnasRNAhybrdScrs[mirna].append(frEnrgy)
					else:
						dMirnasRNAhybrdScrs[mirna] = [frEnrgy]
				else:
					continue
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaRNAhCnts = zeros(cntMirnaNms,dtype = float32)
		aMirnaRNAhEngy = zeros(cntMirnaNms,dtype = float32)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirnasRNAhybrdScrs.has_key(mirna):
				lMirnas = dMirnasRNAhybrdScrs.pop(mirna)
				aMirnaRNAhCnts[mirPos] = len(lMirnas)
				aMirnaRNAhEngy[mirPos] = min(lMirnas)
			else:
				aMirnaRNAhEngy[mirPos] = nan
		return aMirnaRNAhCnts,aMirnaRNAhEngy
	#----------------------------
	#Parse RNAhybrid file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		RNAhybrdRsltsFl = os.path.join(RNAhybrdRsltsFldr, \
		''.join([dt.name,extnsnRNAhybrd]))
		if os.path.exists(RNAhybrdRsltsFl):
			aMirnaRNAhCnts,aMirnaRNAhEngy = prsRNAhybrdFile \
			(RNAhybrdRsltsFl,srtdMirnaFms)
		else:
			aMirnaRNAhCnts = zeros(cntMirnaNms,dtype = float32)
			aMirnaRNAhEngy = zeros(cntMirnaNms,dtype = float32)
			aMirnaRNAhEngy.fill(nan)
		dt.aMirnaRNAhCnts = aMirnaRNAhCnts
		dt.aMirnaRNAhEngy = aMirnaRNAhEngy
	return 0


########################################################
#~ Parse SVMicro results
########################################################
def prsSVMicroRslts(SVMicroRsltsFldr,srtdMirnaFms,extnsnSVMicro,lDtObjct):
	"""
	Input: SVMicroRsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted list of 
	miRNAs of interest. extnsnSVMicro is the extension of the files in 
	the results folder. lDtObjct is a list of data objects, must be 
	the same as that of the files in SVMicroRsltsFldr excepting for the 
	extension. The objects in this list must include the pointer "pos", 
	"name" and "aMirnaNms" (in a congruent way with the input list). 
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", "aMirnaSVMicroCnts", "aMirnaSVMicroScrs"
	(in the same order as "aMirnaNms") for the genes/lncrnas of interest 
	in SVMicroRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsSVMicroFile(SVMicroRsltsFl,srtdMirnaFms):
		"""
		Input: "SVMicroRsltsFl" is an output file from SVMicro program. 
		"srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaSVMicroCnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". aMirnaSVMicroScrs are scores of non-overlapping 
		MREs for miRNAs in "srtdMirnaFms".
		"""	
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaSVMicroCnts = zeros(cntMirnaNms,dtype = float32)
		aMirnaSVMicroScrs = zeros(cntMirnaNms,dtype = float32)	
		#
		dMirnaPos = dict([(mirna,pos) for pos,mirna in \
		enumerate(srtdMirnaFms)])
		lBLckMirnInf = [inf.splitlines() for inf in \
		open(SVMicroRsltsFl).read().split('>') if inf.strip()]
		for eMirna in lBLckMirnInf:
			lPBST = []#info of Potential Binding Site Table
			UTRscr = nan#UTR score
			PBST = False
			hdr = True
			mirna = False
			for eline in eMirna:
				if eline.strip():
					if hdr:
						mirna = eline.split()[0]
						hdr = False
					elif PBST:
						if eline.find \
							('- UTR Feature Table -') > -1:
							PBST = False
						else:
							lPBST.append(eline)
					else:
						if eline.find \
						('- Potential Binding Site Table -') > -1:
							PBST = True
						elif eline.find('UTR score:') > -1:
							UTRscr = float32(eline.split(':')[1]. \
							strip())
			cnts = len(lPBST)-1
			mirPos = dMirnaPos[mirna]
			if cnts > 0:
				aMirnaSVMicroCnts[mirPos] = cnts
			if UTRscr > 0:#only scores > 0 are targets
				aMirnaSVMicroScrs[mirPos] = UTRscr
			else:
				aMirnaSVMicroScrs[mirPos] = nan
		return aMirnaSVMicroCnts,aMirnaSVMicroScrs
	#----------------------------
	#Parse SVMicro file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		SVMicroRsltsFl = os.path.join(SVMicroRsltsFldr, \
		''.join([dt.name,extnsnSVMicro]))
		if os.path.exists(SVMicroRsltsFl):
			aMirnaSVMicroCnts,aMirnaSVMicroScrs = \
			prsSVMicroFile(SVMicroRsltsFl,srtdMirnaFms)
		else:
			aMirnaSVMicroCnts = zeros(cntMirnaNms,dtype = float32)
			aMirnaSVMicroScrs = zeros(cntMirnaNms,dtype = float32)
			aMirnaSVMicroScrs.fill(nan)
		dt.aMirnaSVMicroCnts = aMirnaSVMicroCnts
		dt.aMirnaSVMicroScrs = aMirnaSVMicroScrs
	return 0


########################################################
#~ Parse TargetMiner results
########################################################
def prsTrgtMnrRslts(trgtMnrRsltsFldr,srtdMirnaFms,extnsnTrgtMnr, \
	lDtObjct):
	"""
	Input: trgtMnrRsltsFldr is a folder with files named as the gene 
	containing the results of interest. srtdMirnaFms is a sorted list of 
	miRNAs of interest. extnsnTrgtMnr is the extension of the files in 
	the results folder. lDtObjct is a list of data objects, must be the 
	same as that of the files in trgtMnrRsltsFldr excepting for the 
	extension. The objects in this list must include the pointer "pos", 
	"name" and "aMirnaNms" (in a congruent way with the input list). 
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms", and "aMirnaTrgtMnrCnts" (in the same 
	order as "aMirnaNms") for the genes/lncrnas of interest in 
	trgtMnrRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsTrgtMnrFile(trgtMnrRsltsFl,srtdMirnaFms):
		"""
		Input: "trgtMnrRsltsFl" is an output file from TargetMiner 
		program. "srtdMirnaFms" is a sorted list of miRNAs of interest.
		Output: aMirnaTrgtMnrCnts are the counts of MREs for miRNAs in 
		"srtdMirnaFms". They are in the same order as in "srtdMirnaFms".
		"""
		sPosCnts = set([2,4,6,8])#positions with counts for each MRE kind
		dMirnasTrgtMnrCnts = {}
		hdr = True
		for eline in open(trgtMnrRsltsFl,'r'):
			if hdr:
				hdr = False
			elif eline.strip():
				allinf = eline.split('\t')
				mirnNm = allinf[0]
				mirna = '-'.join(mirnNm.split('-')[1:])
				cntMtchs = sum([int(allinf[p]) for p in sPosCnts if \
				allinf[p]!='-'])
				dMirnasTrgtMnrCnts[mirna] = cntMtchs
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaTrgtMnrCnts = zeros(cntMirnaNms,dtype = float32)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirnasTrgtMnrCnts.has_key(mirna):
				cntMtchs = dMirnasTrgtMnrCnts.pop(mirna)
				aMirnaTrgtMnrCnts[mirPos] = cntMtchs
		return aMirnaTrgtMnrCnts
	#----------------------------
	#Parse TargetMiner file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		trgtMnrRsltsFl = os.path.join(trgtMnrRsltsFldr, \
		''.join([dt.name,extnsnTrgtMnr]))
		if os.path.exists(trgtMnrRsltsFl):
			aMirnaTrgtMnrCnts =	prsTrgtMnrFile(trgtMnrRsltsFl, \
			srtdMirnaFms)
		else:
			aMirnaTrgtMnrCnts = zeros(cntMirnaNms,dtype = float32)
		dt.aMirnaTrgtMnrCnts = aMirnaTrgtMnrCnts
	return 0


########################################################
#~ Parse TargetScan results
########################################################
def prsTrgtScnRslts(trgtScnRsltsFldr,sppMltzCd,srtdMirnaFms, \
	extnsnTrgtScn,lDtObjct):
	"""
	Input: trgtScnRsltsFldr is a folder with files named as the gene 
	containing the results of interest. sppMltzCd is the species code of 
	interest for MultiZ alignment (is a number). srtdMirnaFms is a 
	sorted list of miRNAs of interest. extnsnTrgtScn is the extension of 
	the files in the results folder. lDtObjct is a list of data 
	objects, must be the same as that of the files in trgtScnRsltsFldr 
	excepting for the extension. The objects in this list must include 
	the pointer "pos", "name" and "aMirnaNms" (in a congruent way with 
	the input list).
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMirnaNms" and "aMirnaCnts" (in the same order as 
	"aMirnaNms") for the genes/lncrnas of interest in trgtScnRsltsFldr.
	NOTE: Objects with no results files are going to return scores/
	counts == nan.
	"""
	def prsTrgtScnFile(trgtScnRsltsFl,sppMltzCd,srtdMirnaFms):
		"""
		Input: trgtScnRsltsFl is an output file from targetscan_70.pl 
		script. sppMltzCd is the species code of interest for MultiZ 
		alignment (is a number). srtdMirnaFms is a sorted list of 
		miRNAs of interest. 
		Output: aMirnaCnts are the counts of non-overlapping MREs for 
		miRNAs in srtdMirnaFms found for the species sppMltzCd. They 
		are in 	the same order as in srtdMirnaFms.
		"""
		dMirnasCnt = {}
		notHdr = False
		for eline in open(trgtScnRsltsFl,'r'):
			if notHdr:#skip header
				if eline.strip():
					allinf = eline.split()
					mirna = allinf[1]
					spp = allinf[2]
					if spp == sppMltzCd:
						if dMirnasCnt.has_key(mirna):
							dMirnasCnt[mirna] += 1
						else:
							dMirnasCnt[mirna] = 1
			else:
				notHdr=True
		cntMirnaNms = len(srtdMirnaFms)
		aMirnaCnts = zeros(cntMirnaNms,dtype = float32)
		for mirPos in xrange(cntMirnaNms):
			mirna = srtdMirnaFms[mirPos]
			if dMirnasCnt.has_key(mirna):
				aMirnaCnts[mirPos] = dMirnasCnt[mirna]
		return aMirnaCnts
	#----------------------------
	#Parse TargetScan file
	cntMirnaNms = len(srtdMirnaFms)
	for dt in lDtObjct:
		trgtScnRsltsFl = os.path.join(trgtScnRsltsFldr,''.join([dt.name, \
		extnsnTrgtScn]))
		if os.path.exists(trgtScnRsltsFl):
			aMirnaCnts = prsTrgtScnFile(trgtScnRsltsFl,sppMltzCd, \
			srtdMirnaFms)
		else:
			aMirnaCnts = zeros(cntMirnaNms,dtype = float32)
		dt.aMirnaCnts = aMirnaCnts	
	return 0


########################################################
#~ Parse ENSEMBL annotations (obtain with rtrnAnntTrscrptsSubSt/
# rtrnAnntGenesSubSt) and make a dictionary of names and annotation
########################################################
def mkDctAnnt(infNmAnnt,clmnNme,clmnAnnt,excldNs=True):
	"""
	Input: infNmAnnt is a file with at least three tab-delimted 
	columns. clmnNme is the column with the name (0-based). clmnAnnt
	is the clmn with the annotation (0-based). Optionally, excldNs
	is a boolean to choose if null anotations 'N' are excluded 
	(default=True).
	Output: dNameAnnt is a dicitonary with the name as key and the 
	annotations as values.
	"""
	#----------------------------
	# Work with file
	dNameAnnt={}
	for l in open(infNmAnnt,'r'):
		if l.strip():
			l=l.splitlines()[0].split('\t')
			name=l[clmnNme]
			annt=l[clmnAnnt].replace(' ','_')
			if excldNs: 
				if annt!='N':
					try:
						dNameAnnt[name].add(annt)
					except:
						dNameAnnt[name]=set([annt])
			else:
				try:
					dNameAnnt[name].add(annt)
				except:
					dNameAnnt[name]=set([annt])
	#----------------------------
	# Merge annotations for the same name
	dNameAnnt=dict([(name,'.'.join(sorted(dNameAnnt[name]))) for name \
	in dNameAnnt.keys()])
	return dNameAnnt

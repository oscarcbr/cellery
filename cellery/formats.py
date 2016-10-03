#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  formats.py part of cellery (ceRNAs linking inference)
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
Methods to work with different formats
"""

########################################################
#~ Import libraries.
########################################################
from itertools import product
from multiprocessing import Queue,Process
from numpy import empty,array

import os

########################################################
#~ Format TargetScan miRNA file
########################################################
def mkArryFlsMirnaTrgtScn(fstaFlmirbaseFrmtd,tgrtScnFlmrnFrmtd, \
	mirbaseSpp,sppMltzCd):
	"""
	Input: fstaFlmirbaseFrmtd is the full path to a Fasta-file of
	miRNAs of interest with header in format of miRBase. 
	tgrtScnFlmrnFrmtd is the full path to the save the miRNA input file 
	in TargetScan format. mirbaseSpp is the species code in miRBase base 
	(i.e. ==KEGG format: {hsa,mmu}). sppMltzCd is the species code of 
	interest for MultiZ alignment (is a number).
	Output: write the TargetScan-formated miRNA file.
	"""
	def frmtFastaTOTrgtScn(fstaFlmirbaseFrmtd,mirbaseSpp,sppMltzCd):
		"""
		"""
		mirnFrmtdFlTgtScn = []
		for eseq in open(fstaFlmirbaseFrmtd).read().split('>'):
			if eseq.strip():
				hdr,seq = eseq.splitlines()
				hdrInf = hdr.split()
				mirnaNm,fmlyID = hdrInf[0],hdrInf[-1]
				spp = mirnaNm.split('-')[0]
				if spp==mirbaseSpp:
					frmtdInf = '\t'.join([fmlyID,seq[1:8],sppMltzCd])
					mirnFrmtdFlTgtScn.append(frmtdInf)
		mirnFrmtdFlTgtScn.sort()
		frtmdMirnaTgtScn='\n'.join(mirnFrmtdFlTgtScn)
		return frtmdMirnaTgtScn
	#----------------------------
	#~ Convert fasta to TargetScan format
	frtmdMirnaTgtScn = frmtFastaTOTrgtScn(fstaFlmirbaseFrmtd,mirbaseSpp, \
	sppMltzCd)
	flmirbaseFrmtd = open(tgrtScnFlmrnFrmtd,'w')
	flmirbaseFrmtd.write(frtmdMirnaTgtScn)
	flmirbaseFrmtd.close()
	return 0


########################################################
#~ Format TargetMiner pairwise file
########################################################
def mkPrwsFlTrgtMnr(lDtObjcts,lMirnasGrps,pthFldrPrwsFrmtdFlsTrgtMnr, \
	extnsnPrwsTrgtMnr,nThrds):
	"""
	Input: lDtObjcts is a list of data objects with pointer "name". 
	lMirnasGrps is the group of mirna headers limited to 999. 
	pthFldrPrwsFrmtdFlsTrgtMnr is the path to write the tabular- 
	formated files to input in TargetMiner. extnsnPrwsTrgtMnr is the 
	extension of the tabular-formated files to input in TargetMiner. 
	nThrds is the number to run in multithread.
	Output: Writes a tabular-formated file with all the possible 
	pairwise comparison between mirnas and genes.
	"""
	def mthrdWrtPrws(qInJobs,qOutRslts,pthFldrPrwsFrmtdFlsTrgtMnr, \
		extnsnPrwsTrgtMnr,lMirnasGrps):
		"""
		"""
		for dtNm,lMirnasPos in iter(qInJobs.get,'STOP'):
			opndPrwsFrmtdFlsTrgtMnr = os.path.join( \
			pthFldrPrwsFrmtdFlsTrgtMnr, ''.join([dtNm,'_%s'%lMirnasPos, \
			extnsnPrwsTrgtMnr]))
			outLst =['\t'.join([mirnaNm.split()[0],dtNm]) for mirnaNm in \
			lMirnasGrps[lMirnasPos]]
			outLst.append('>')
			opndPrwsFrmtdFlsTrgtMnr = open(opndPrwsFrmtdFlsTrgtMnr,'w')
			opndPrwsFrmtdFlsTrgtMnr.write('\n'.join(outLst))
			opndPrwsFrmtdFlsTrgtMnr.close()
			qOutRslts.put(0)
	#
	qInJobs = Queue()
	cntMirnasGrps = len(lMirnasGrps)
	for dt,lMirnasPos in product(lDtObjcts,xrange(cntMirnasGrps)):
		dtNm = dt.name
		qInJobs.put((dtNm,lMirnasPos))
	#write in multithread
	totalPrs = cntMirnasGrps*len(lDtObjcts)
	qOutRslts = Queue()
	for p in xrange(nThrds):
		Process(target=mthrdWrtPrws,args=(qInJobs,qOutRslts, \
		pthFldrPrwsFrmtdFlsTrgtMnr,extnsnPrwsTrgtMnr,lMirnasGrps)).start()
	for p in xrange(totalPrs):
		qOutRslts.get()
	for p in xrange(nThrds):
		qInJobs.put('STOP')
	return 0


########################################################
#~ Format TargetMiner miRNA file
########################################################
def procssMirnTrgtMnr(pthFldrFstaMirnsFls,fstasMirnaBsFl, \
	mirnaFrmtdFlTrgtMnr,pthFldrMirnsFrmtdFlsTrgtMnr,extnsnMirnaTrgtMnr):
	"""
	Input: fstasMirnaBsFl is the full path to the fasta-formated mirna 
	sequences. pthFldrFstaMirnsFls is the full path to fstasMirnaBsFl.  
	mirnaFrmtdFlTrgtMnr is the output file to write the TargetMiner-
	formated output file. pthFldrMirnsFrmtdFlsTrgtMnr is the full path
	to mirnaFrmtdFlTrgtMnr. extnsnMirnaTrgtMnr is the extension of the 
	TargetMiner-formated output file.
	Output: Files mirnaFrmtdFlTrgtMnr_cnt with mirna sequences in 
	TargetMiner format, limiting the number to 999. lMirnasGrpsOut is 
	the group of mirna headers limited to 999.
	"""
	fstasMirnaBsFlNm = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	lMirnasGrps = []
	lMirnas = []
	cnt = -1
	for eline in open(fstasMirnaBsFlNm,'r'):
		cnt += 1
		if cnt == 1000:
			lMirnas.sort()
			lMirnasGrps.append(lMirnas)
			lMirnas = []
			cnt = 0
		if eline.strip():
			if eline[0] == '>':
				hdr = ' '.join(eline.split()[:2])
			else:
				seq = eline.strip()
				lMirnas.append((hdr,seq))
	#
	if cnt<1000:
		lMirnas.sort()
		lMirnasGrps.append(lMirnas)
	#
	lMirnasGrpsOut = []
	cnt = -1
	for lMirnas in lMirnasGrps:
		cnt += 1
		mirnaFrmtdFlTrgtMnrNm = os.path.join(pthFldrMirnsFrmtdFlsTrgtMnr, \
		''.join([mirnaFrmtdFlTrgtMnr,'_%s'%cnt,extnsnMirnaTrgtMnr]))
		opndMirnaFrmtdFlTrgtMnrNm = open(mirnaFrmtdFlTrgtMnrNm,'w')
		lMirnasOut = []
		for hdr,seq in lMirnas:
			lMirnasOut.append(hdr)
			opndMirnaFrmtdFlTrgtMnrNm.write('\n'.join([hdr,seq,'']))
		opndMirnaFrmtdFlTrgtMnrNm.write('>')
		opndMirnaFrmtdFlTrgtMnrNm.close()
		lMirnasGrpsOut.append(lMirnasOut)
	#
	return lMirnasGrpsOut


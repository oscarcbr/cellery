#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  wrapper.py part of cellery (ceRNAs linking inference)
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
Libraries with wrappers to run MRE predicting algorithms for cellery
"""


########################################################
#~ Import external libraries
########################################################
from cellery import alignments
from cellery import exceptions
from cellery import formats
from cellery import parsers
from itertools import product
from multiprocessing import Queue,Process
from numpy import array
from tempfile import NamedTemporaryFile

import os,shutil,BeautifulSoup

########################################################
#~ Add predictions to a list of objects from another list of objects
########################################################
def addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd,lAttrbtsToAdd):
	"""
	Input: lObjcts_trgt is the list of objects that will receive the 
	attributes in the list lAttrbtsToAdd from lObjcts_toAdd. 
	lObjcts_toAdd is the list of objects that will donate the attributes 
	in the list lAttrbtsToAdd to lObjcts_trgt. lAttrbtsToAdd is a list 
	of attributes to be transferred. 
	Output: lObjcts_trgt will contain the objects with attributes in the 
	list lAttrbtsToAdd donated from lObjcts_toAdd updated.
	"""
	assert len(lObjcts_trgt)==len(lObjcts_toAdd)
	lenObjcts = len(lObjcts_trgt)
	for pos in xrange(lenObjcts):
		assert lObjcts_trgt[pos].name==lObjcts_toAdd[pos].name
		for atrbt in lAttrbtsToAdd:
			atrbtVl = getattr(lObjcts_toAdd[pos],atrbt)
			setattr(lObjcts_trgt[pos],atrbt,atrbtVl)


########################################################
#~ Run mirMap in multithread
########################################################
def queueMirMap(aFstfrmtdFls,mrnMirmapFrmtdFl,pthFldrUTRFrmtdFls, \
	pthFldrmrnMirmapFrmtdFls,pthOutFldrMirMap,pthTmpFldr,pthPython, \
	pthPythonLib,pthPhast,pthMirMapLib,extnsnMirMap,extnsnUTRMirMap, \
	phastModFl,phastTree,nPrlProcss,nThrds):
	"""
	Input: aFstfrmtdFls is an array of UTR Fasta-formated file names. 
	mrnMirmapFrmtdFl is an array of miRNA Fasta-formated file names. 
	pthFldrUTRFrmtdFls is a path to the folder with the files in 
	aUTRFrmtdFls. pthFldrmrnMirmapFrmtdFls is a path to the folder with 
	the files in mrnMirmapFrmtdFl. pthOutFldrMirMap is a path to the 
	folder to save the MirMap results. pthTmpFldr is a path to a 
	temporal folder to save the queue input files. pthPython is the full
	path to python. pthPythonLib is the full path to the python library.
	pthPhast is the full path to the phast executable. pthMirMapLib is 
	the full path to the mirMap python library. extnsnMirMap is the 
	extension of mirMap file results. extnsnUTRMirMap is the extension 
	for the UTRs of the objects. phastModFl is the phastCons model file. 
	phastTree is the phastCons tree file (can be obtained from the UCSC 
	browser multiZ version of interest). nPrlProcss is the number of 
	parallel processes to run mirMap. nThrds is the number of cores to 
	run.
	Output: The results of MirMap written in the pthOutFldrMirMap. 
	sOutFlsMirMap is a the set of names for these files.
	NOTE: in the input file the reference sequence MUST be the first 
	one, as default in alignments produced by cellery.
	"""
	#
	def rtrnMirmapScrpt(UTRseq,mrnMirmapFrmtdFlNm,UTRFrmtdFlNm,phastTree, \
		phastModFl,pthMirMapLib,pthPhast,outFl):
		"""
		return mirmap script with default parameters
		"""
		scrpt = '\n'.join([
		'import mirmap,os',
		'import mirmap.library_link',
		'#from numpy import array',
		'seq_target = "%s"'%UTRseq,
		'aMirnsSeqs = []',
		'aMirnsNames = []',
		'for echl in open("%s","r"):'%mrnMirmapFrmtdFlNm,
		'\tif echl.strip():',
		'\t\tif echl[0] == ">":',
		'\t\t\taMirnsNames.append("_".join(echl[1:].split()[0].split("-")[1:]))',
		'\t\telse:',
		'\t\t\taMirnsSeqs.append(echl.strip())',
		'outL=[]',
		'for mrnMirmapFrmtdFlpos in xrange(len(aMirnsSeqs)):',
		'\tmirnaSeq = aMirnsSeqs[mrnMirmapFrmtdFlpos]',
		'\tseq_mirna = "%s"%mirnaSeq',
		'\tmim = mirmap.mm(seq_target, seq_mirna)',
		'\tmim.find_potential_targets_with_seed()',
		'\tmim.end_sites',
		'\tif mim.end_sites:',
		'\t\tmrnMirmapFrmtdFlPrfx = aMirnsNames[mrnMirmapFrmtdFlpos]',
		'\t\tmim.eval_tgs_au()',
		'\t\tmim.eval_tgs_pairing3p()',
		'\t\tmim.eval_tgs_position()',
		'\t\tmim.prob_binomial',
		'\t\tmim.libs = mirmap.library_link.LibraryLink("%s")'% \
		pthMirMapLib,
		'\t\tmim.exe_path = mirmap.library_link.LibraryLink("%s")'% \
		pthMirMapLib,
		'\t\tmim.dg_duplex',
		'\t\tmim.dg_open',
		'\t\tmim.prob_exact',
		'\t\tlib = mirmap.if_clib_phast.Phast(library_path="%s")'% \
		pthMirMapLib,
		'\t\tmim.eval_selec_phylop(aln_fname="%s",pathphast="%s",mod_fname="%s")'% \
		(UTRFrmtdFlNm,pthPhast,phastModFl),
		'\t\tmim.eval_cons_bls(aln_fname="%s",pathphast="%s",subst_model="REV",tree="%s")'% \
		(UTRFrmtdFlNm,pthPhast,phastTree),
		'\t\tmim.eval_score()',
		'\t\tmirmapscr = mim.score',
		'\t\toutRep = mim.report().__str__()',
		'\t\tif mirmapscr:',
		'\t\t\toutL.append(">%s\\n%s\\nmirMap_Score:\t%s"%(mrnMirmapFrmtdFlPrfx,outRep,mirmapscr))',
		'\t\telse:',
		'\t\t\toutL.append(">%s\\n%s"%(mrnMirmapFrmtdFlPrfx,outRep))',
		'#',
		'outFile = open("%s","w")'%outFl,
		"outFile.write('\\n\\n'.join(outL))",
		'outFile.close()'])
		return scrpt
	#
	def wrtTmp(tmpFl,pthPython,UTRseq,mrnMirmapFrmtdFlNm,UTRFrmtdFlNm, \
		phastTree,phastModFl,pthMirMapLib,pthPhast,pthPythonLib,outFl):
		"""
		Write a temporary file to run mirMap in a queue
		"""
		#write mirmap script
		fScrpt = '%s.py'%os.path.splitext(tmpFl)[0]
		fOpnScrpt = open(fScrpt,'w')
		scprt = rtrnMirmapScrpt(UTRseq,mrnMirmapFrmtdFlNm,UTRFrmtdFlNm, \
		phastTree,phastModFl,pthMirMapLib,pthPhast,outFl)
		fOpnScrpt.write(scprt)
		fOpnScrpt.close()
		#write tmp bash file
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nexport PATH=%s:%s:$PATH\n'% \
		(pthPython,pthPhast))
		fOpn.write('export PYTHONPATH=%s\n'%pthPythonLib)
		fOpn.write('python %s > /dev/null'%fScrpt)
		fOpn.close()
	#
	def mthrdWrtTmp(qInJobs,qOutRslts,mrnMirmapFrmtdFlNm, \
		pthPython,phastTree,phastModFl,pthMirMapLib,pthPhast, \
		pthPythonLib):
		"""
		Multithread writing a temporary file to run mirMap in a queue
		"""
		for UTRFrmtdFlNm,outFl,tmpFl in iter(qInJobs.get, \
			'STOP'):
			UTRseq = ''.join([l for l in open(UTRFrmtdFlNm).read() \
			.split('>')[1].splitlines()[1:] if l.strip()]).replace('-',\
			'')
			wrtTmp(tmpFl,pthPython,UTRseq,mrnMirmapFrmtdFlNm, \
			UTRFrmtdFlNm,phastTree,phastModFl,pthMirMapLib,pthPhast, \
			pthPythonLib,outFl)
			qOutRslts.put(0)
	#
	assert pthPython and pthPythonLib and pthPhast and pthMirMapLib
	#
	mrnMirmapFrmtdFlNm = os.path.join(pthFldrmrnMirmapFrmtdFls, \
	mrnMirmapFrmtdFl)
	sOutFlsMirMap = set()		
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	qInJobs = Queue()
	#
	for fstFrmtdFl in aFstfrmtdFls:
		UTRFrmtdFlPrfx = os.path.splitext(fstFrmtdFl)[0]
		UTRFrmtdFl = ''.join([UTRFrmtdFlPrfx,extnsnUTRMirMap])
		outFlNm = ''.join([UTRFrmtdFlPrfx,extnsnMirMap])
		outFl = os.path.join(pthOutFldrMirMap,outFlNm)
		sOutFlsMirMap.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			UTRFrmtdFlNm = os.path.join(pthFldrUTRFrmtdFls,UTRFrmtdFl)
			qInJobs.put((UTRFrmtdFlNm,outFl,tmpFl))
	#write in multithread
	qOutRslts = Queue()
	for p in xrange(nThrds):
		Process(target=mthrdWrtTmp,args=(qInJobs,qOutRslts, \
		mrnMirmapFrmtdFlNm,pthPython,phastTree,phastModFl,pthMirMapLib, \
		pthPhast,pthPythonLib)).start()
	cnt=0
	for p in xrange(cntTmp):
		qOutRslts.get()
	for p in xrange(nThrds):
		qInJobs.put('STOP')
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in os.listdir \
		(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsMirMap


########################################################
#~ Run Miranda in multithread
########################################################
def queueRunMiranda(aFstfrmtdFls,mrnMirandaFrmtdFl,pthFldrFstfrmtdFls, \
	pthFldrMirnsFrmtdFls,pthOutFldrMiranda,pthTmpFldr,pthMiranda, \
	extnsnMiranda,nPrlProcss):
	"""
	Input: aFstfrmtdFls is an array of UTR Fasta-formated file names. 
	mrnMirandaFrmtdFl is an array of miRNA Fasta-formated file name. 
	pthFldrFstfrmtdFls is a path to the folder with the files in 
	aFstfrmtdFls. pthFldrMirnsFrmtdFls is a path to the folder with the 
	files in mrnMirandaFrmtdFl. pthOutFldrMiranda is a path tothe folder 
	to save the miRanda results. pthTmpFldr is a path to a temporal 
	folder to save the queue input files. pthMiranda is the path to the
	miranda executable. extnsnMiranda is the extension of miRanda file 
	results. nPrlProcss is the number of parallel processes to run 
	miRanda.
	Output: The results of miRanda written in the pthOutFldrMiranda. 
	sOutFlsMiranda is a the set of names for these files.
	"""
	#
	def wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthMiranda):
		"""
		Write a temporary file to run Miranda in a queue
		"""
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nPATH=%s:$PATH\n'%pthMiranda)
		fOpn.write('miranda %s %s -quiet -out %s'%(pthMirnsFrmtdFl, \
		pthFstFrmtdFl,outFl))
		fOpn.close()
	#
	assert pthMiranda
	#
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	sOutFlsMiranda = set()
	cntTmp = 0
	pthMirnsFrmtdFl = os.path.join(pthFldrMirnsFrmtdFls, \
	mrnMirandaFrmtdFl)
	for fstFrmtdFl in aFstfrmtdFls:
		fstFrmtdFlPrfx = os.path.splitext(fstFrmtdFl)[0]
		outFlNm = ''.join([fstFrmtdFlPrfx,extnsnMiranda])
		outFl = os.path.join(pthOutFldrMiranda,outFlNm)
		sOutFlsMiranda.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			pthFstFrmtdFl = os.path.join(pthFldrFstfrmtdFls,fstFrmtdFl)
			wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthMiranda)
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in os.listdir(pthTmpFldr) \
		if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsMiranda


########################################################
#~ Run PITA in multithread
########################################################
def queueRunPITA(aUTRFrmtdFls,mrnPITAFrmtdFl,pthFldrfstFrmtdFls, \
	pthFldrMirnsFrmtdFls,pthOutFldrPITA,pthTmpFldr,pthPITA,extnsnPITA, \
	nPrlProcss):
	"""
	Input: aUTRFrmtdFls is an array of UTR PITA-formated file names. 
	mrnPITAFrmtdFl is a miRNA fasta-formated file names. 
	pthFldrFstfrmtdFls is a path to the folder with the files in 
	aFstfrmtdFls. pthFldrMirnsFrmtdFls is a path to the folder with the 
	files in mrnPITAFrmtdFl. pthOutFldrPITA is a path to the folder to 
	save the PITA results. pthTmpFldr is a path to a temporal folder to 
	save the queue input files. pthPITA is the path the PITA executable.
	extnsnPITA is the extension of PITA file results. nPrlProcss is the 
	number of parallel processes to run PITA.
	Output: The results of PITA written in the pthOutFldrPITA. 
	sOutFlsPITA is a the set of names for these files.
	"""
	#
	def wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthPITA):
		"""
		Write a temporary file to run PITA in a queue
		"""
		brtOutFl = os.path.splitext(outFl)[0]
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nPATH=%s:$PATH\n'%pthPITA)
		fOpn.write('pita_prediction.pl -utr %s -mir %s -prefix %s\n' \
		%(pthFstFrmtdFl,pthMirnsFrmtdFl,brtOutFl))
		fOpn.write('mv %s_pita_results_targets.tab %s\n'%(brtOutFl, \
		outFl))
		fOpn.write('rm %s_*\n'%brtOutFl)
		fOpn.close()
	#
	assert pthPITA
	#
	sOutFlsPITA = set()
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	pthMirnsFrmtdFl = os.path.join(pthFldrMirnsFrmtdFls, \
	mrnPITAFrmtdFl)
	for fstFrmtdFl in aUTRFrmtdFls:
		fstFrmtdFlPrfx = os.path.splitext(fstFrmtdFl)[0]
		outFlNm = ''.join([fstFrmtdFlPrfx,extnsnPITA])
		outFl = os.path.join(pthOutFldrPITA,outFlNm)
		sOutFlsPITA.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			pthFstFrmtdFl = os.path.join(pthFldrfstFrmtdFls,fstFrmtdFl)
			wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthPITA)
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in os.listdir \
		(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsPITA


########################################################
#~ Run RNAhybrid in multithread
########################################################
def queueRunRNAhybrd(aFstfrmtdFls,mrnRNAhybrdFrmtdFl, \
	pthFldrFstfrmtdFls,pthFldrMirnsFrmtdFls,pthOutFldrRNAhybrd, \
	pthTmpFldr,pthRNAhybrd,extnsnRNAhybrd,nPrlProcss):
	"""
	Input: aFstfrmtdFls is an array of UTR Fasta-formated file names. 
	mrnRNAhybrdFrmtdFl is a miRNA fasta-formated file. 
	pthFldrFstfrmtdFls is a path to the folder with the files in 
	aFstfrmtdFls. pthFldrMirnsFrmtdFls is a path to the folder with the 
	files in mrnRNAhybrdFrmtdFl. pthOutFldrRNAhybrd is a path to the 
	folder to save the RNAhybrid results. pthTmpFldr is a path to a 
	temporal folder to save the queue input files. pthRNAhybrd is the 
	path to the RNAhybrid executable. extnsnRNAhybrd is the extension of 
	RNAhybrid file results. nPrlProcss is the number of parallel 
	processes to run RNAhybrid.
	Output: The results of RNAhybrid written in the pthOutFldrRNAhybrd. 
	sOutFlsRNAhybrd is a the set of names for these files.
	"""
	#	
	def wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthRNAhybrd, \
		seqLen,minEnrgy='-20'):
		"""
		Write a temporary file to run RNAhybrid in a queue
		"""
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nPATH=%s:$PATH\n'%pthRNAhybrd)
		fOpn.write('RNAhybrid -f 2,8 -s 3utr_human -t %s -q %s -e %s -c \
		-m %s > %s'%(pthFstFrmtdFl,pthMirnsFrmtdFl,minEnrgy,seqLen, \
		outFl))
		fOpn.close()
	#
	assert os.path.exists(pthRNAhybrd)
	#
	sOutFlsRNAhybrd = set()
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	pthMirnsFrmtdFl = os.path.join(pthFldrMirnsFrmtdFls, \
	mrnRNAhybrdFrmtdFl)
	for fstFrmtdFl in aFstfrmtdFls:
		fstFrmtdFlPrfx = os.path.splitext(fstFrmtdFl)[0]
		outFlNm = ''.join([fstFrmtdFlPrfx,extnsnRNAhybrd])
		outFl = os.path.join(pthOutFldrRNAhybrd,outFlNm)
		sOutFlsRNAhybrd.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			pthFstFrmtdFl = os.path.join(pthFldrFstfrmtdFls,fstFrmtdFl)
			seqLen = len(''.join([l for l in open(pthFstFrmtdFl).read() \
			.splitlines()[1:] if l.strip()]))
			wrtTmp(pthFstFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl, \
			pthRNAhybrd,seqLen)
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in \
		os.listdir(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsRNAhybrd


########################################################
#~ Run SVMicro in multithread
########################################################
def queueRunSVMicro(aUTRFrmtdFls,mrnSVMicroFrmtdFl,pthFldrUTRFrmtdFls, \
	pthFldrmrnSVMicroFrmtdFls,pthOutFldrSVMicro,pthTmpFldr,pthSVMicro, \
	extnsnSVMicro,nPrlProcss,nThrds):
	"""
	Input: aUTRFrmtdFls is an array of UTR SVMicro-formated file names. 
	mrnSVMicroFrmtdFl is a miRNA fasta-formated file names. 
	pthFldrUTRFrmtdFls is a path to the folder with the UTR SVMicro-
	formated files. pthFldrmrnSVMicroFrmtdFls is a path to the folder 
	with the files in mrnSVMicroFrmtdFl. pthOutFldrSVMicro is a path to 
	the folder with the SVMicro results. pthTmpFldr is a path to a 
	temporal folder to save the queue input files. pthSVMicro is the 
	path to the SVMicro executable. extnsnSVMicro is the extension of 
	SVMicro file results. nPrlProcss is the number of parallel processes 
	to run SVMicro. nThrds is the number of threads in multithreading. 
	Output: The results are SVMicro results written in the 
	pthOutFldrSVMicro. sOutFlsSVMicro is a the set of names for these 
	files.
	"""
	#
	def rtrnSVNMicroScrpt(UTRFrmtdFlPrfx,mrnSVMicroFrmtdFlNm,tmpFlScrpt, \
		pthSVMicro,outFl):
		"""
		return SVNMicro script with default parameters
		"""
		scrpt = '\n'.join([
		'import subprocess',
		'aMirnsSeqs = []',
		'aMirnsNames = []',
		'for echl in open("%s","r"):'%mrnSVMicroFrmtdFlNm,
		'\tif echl.strip():',
		'\t\tif echl[0] == ">":',
		'\t\t\taMirnsNames.append("_".join(echl[1:].split()[0].split("-")[1:]))',
		'\t\telse:',
		'\t\t\taMirnsSeqs.append(echl.strip())',
		'#',
		'outL=[]',
		'for mrnSVMicroFrmtdFlpos in xrange(len(aMirnsSeqs)):',
		'\tmirnaSeq = aMirnsSeqs[mrnSVMicroFrmtdFlpos]',
		'\trnSVMicro = subprocess.Popen(["perl","%s","-m",mirnaSeq,"-u","%s"],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd="%s")'% \
		(tmpFlScrpt,UTRFrmtdFlPrfx,pthSVMicro),
		'\toutf = rnSVMicro.communicate()[0]',
		'\tif outf.find("\\nUTR score: -2\\n")==-1 and outf.find("UTR score:")>-1:',
		'\t\tmrnSVMicroFrmtdFlPrfx = aMirnsNames[mrnSVMicroFrmtdFlpos]',
		'\t\toutL.append(">%s\\n%s"%(mrnSVMicroFrmtdFlPrfx,outf))',
		'#',
		'outFile = open("%s","w")'%outFl,
		"outFile.write('\\n\\n'.join(outL))",
		'outFile.close()'])
		return scrpt
	#
	def wrtTmp(UTRFrmtdFlPrfx,pthUTRFrmtdFls,mrnSVMicroFrmtdFlNm,outFl, \
		tmpFl,pthSVMicro,opnSvmicroScrpt):
		"""
		Write a temporary file to run SVMicro in a queue
		"""
		#
		chngdScrpt = ''.join([opnSvmicroScrpt[0],'$utr_consv_file = "', \
		pthUTRFrmtdFls,'";#',opnSvmicroScrpt[1]])
		tmpFlScrpt = ''.join([tmpFl,'.pl'])
		fOpn = open(tmpFlScrpt,'w')
		fOpn.write(chngdScrpt)
		fOpn.close()
		#
		pthnScrpt = rtrnSVNMicroScrpt(UTRFrmtdFlPrfx, \
		mrnSVMicroFrmtdFlNm,tmpFlScrpt,pthSVMicro,outFl)
		tmpFlScrpt = ''.join([tmpFl,'.py'])
		fOpn = open(tmpFlScrpt,'w')
		fOpn.write(pthnScrpt)
		fOpn.close()
		#
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nPATH=%s:$PATH\n'%pthSVMicro)
		fOpn.write('python %s'%tmpFlScrpt)
		fOpn.close()
	#
	assert pthSVMicro
	opnSvmicroScrpt = open(os.path.join(pthSVMicro,'svmicro.pl')).read()
	#correct bugs for long sequences/targets
	opnSvmicroScrpt = opnSvmicroScrpt.replace \
	('`echo -en "$all_site_vectors\\n" >> $svm_sample_file`;', \
	'open my $out,">",$svm_sample_file;\nprint $out %s'% \
	'"$all_site_vectors\\n";\nclose $out;')
	opnSvmicroScrpt = opnSvmicroScrpt.replace \
	('`echo -en "$utr_vector\\n" > $svm_sample_file`;', \
	'open my $outa,">",$svm_sample_file;\nprint $outa "$utr_vector\\n"%s'% \
	';\nclose $outa;')
	opnSvmicroScrpt = opnSvmicroScrpt.split('$utr_consv_file = ')
	#
	mrnSVMicroFrmtdFlNm = os.path.join(pthFldrmrnSVMicroFrmtdFls, \
	mrnSVMicroFrmtdFl)
	sOutFlsSVMicro = set()
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	for UTRFrmtdFl in aUTRFrmtdFls:
		pthUTRFrmtdFls = os.path.join(pthFldrUTRFrmtdFls,UTRFrmtdFl)
		UTRFrmtdFlPrfx = os.path.splitext(UTRFrmtdFl)[0]
		outFlNm = ''.join([UTRFrmtdFlPrfx,extnsnSVMicro])
		outFl = os.path.join(pthOutFldrSVMicro,outFlNm)
		sOutFlsSVMicro.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			wrtTmp(UTRFrmtdFlPrfx,pthUTRFrmtdFls,mrnSVMicroFrmtdFlNm, \
			outFl,tmpFl,pthSVMicro,opnSvmicroScrpt)
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'%(pthTmpFldr, \
		tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in \
		os.listdir(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)	
		return sOutFlsSVMicro


########################################################
#~ Run TargetMiner in multithread
########################################################
def queueRunTrgtMnr(aDtTrgtMnr,cntMirnasGrps,pthFldrUTRFrmtdFlsTrgtMnr, \
	extnsnUTRTrgtMnr,mirnaFrmtdFlTrgtMnr,pthFldrMirnsFrmtdFlsTrgtMnr, \
	extnsnMirnaTrgtMnr,pthFldrPrwsFrmtdFlsTrgtMnr,extnsnPrwsTrgtMnr, \
	pthOutFldrTrgtMnr,pthTmpFldr,pthTrgtMnr,extnsnTrgtMnr,nPrlProcss, \
	nThrds):
	"""
	Input: aDtTrgtMnr is an array of gene/lncrna names of interest. 
	pthFldrUTRFrmtdFlsTrgtMnr is the full path to the UTR files.  
	extnsnUTRTrgtMnr is the extension of the UTR files. cntMirnasGrps is 
	the number of TargetMiner- formated mirna files. mirnaFrmtdFlTrgtMnr 
	is the prefix of the TargetMiner-formated mirna file. 
	pthFldrMirnsFrmtdFlsTrgtMnr is a path to the folder with the files 
	with prefix mirnaFrmtdFlTrgtMnr. extnsnMirnaTrgtMnr is the extension 
	of the TargetMiner-formated mirna file. pthFldrPrwsFrmtdFlsTrgtMnr 
	is a path to the folder with the files in TargetMiner-formated pair 
	files. extnsnPrwsTrgtMnr is the extension of the TargetMiner-formated 
	pair files. pthOutFldrTrgtMnr is the output folder. pthTmpFldr is a 
	path to a temporal folder to temporarly process files. pthTrgtMnr is
	the to the TargetMiner executable. extnsnTrgtMnr is the extension to 
	the output files from TargetMiner. nPrlProcss is the number of 
	parallel processes to run TargetMiner.
	Output: The merged results of TargetMiner for each gene are written 
	in pthOutFldrTrgtMnr. sOutFlsTrgtMnr is a the set of names for these 
	files.
	"""
	#
	def wrtTmp(dtNm,UTRFrmtdFlTrgtMnrNm,cntMirnasGrps, \
		mirnaFrmtdFlTrgtMnr,pthFldrMirnsFrmtdFlsTrgtMnr, \
		extnsnMirnaTrgtMnr,pthFldrPrwsFrmtdFlsTrgtMnr,extnsnPrwsTrgtMnr, \
		outFl,pthTmpFldr,tmpFl,pthTrgtMnr):
		"""
		Write a temporary file to run TargetMiner in a queue
		"""
		tmpFldr = os.path.join(pthTmpFldr,'%s.d'%dtNm)
		trgMnrLiOrg = os.path.join(pthTrgtMnr,'libsvm-2.88')
		trgMnrLiDtny = os.path.join(tmpFldr,'libsvm-2.88')
		trgMnrSVMdlOrgn = os.path.join(pthTrgtMnr,'svm.model')
		trgMnrSVMdlDtny = os.path.join(tmpFldr,'svm.model')
		trgMnrRngOrgn = os.path.join(pthTrgtMnr,'range')
		trgMnrRngDtny = os.path.join(tmpFldr,'range')
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nPATH=%s:$PATH\n'%pthTrgtMnr)
		fOpn.write('mkdir %s\ncd %s\n'%(tmpFldr,tmpFldr))
		fOpn.write('ln -sf %s %s\n'%(trgMnrLiOrg,trgMnrLiDtny))
		fOpn.write('ln -sf %s %s\n'%(trgMnrSVMdlOrgn,trgMnrSVMdlDtny))
		fOpn.write('ln -sf %s %s\n'%(trgMnrRngOrgn,trgMnrRngDtny))
		for cnt in xrange(cntMirnasGrps):
			prwfl = os.path.join(pthFldrPrwsFrmtdFlsTrgtMnr, \
			'%s_%s%s'%(dtNm,cnt,extnsnPrwsTrgtMnr))
			mirnafl = os.path.join(pthFldrMirnsFrmtdFlsTrgtMnr, \
			'%s_%s%s'%(mirnaFrmtdFlTrgtMnr,cnt,extnsnMirnaTrgtMnr))
			fOpn.write('TargetMiner -test %s -3utr %s -mir %s\n' \
			%(prwfl,UTRFrmtdFlTrgtMnrNm,mirnafl))
			orgnlFl = os.path.join(tmpFldr, \
			'TargetMiner_Prediction.html')
			cntFl = os.path.join(tmpFldr,'%s.html'%cnt)
			fOpn.write('mv %s %s\n'%(orgnlFl,cntFl))
		fOpn.close()
	#
	def mrgRslts(dtNm,pthTmpFldr):
		"""
		Merge results from inidividual TargetMiner html result files
		"""
		tmpFldr = os.path.join(pthTmpFldr,'%s.d'%dtNm)
		lRslts = [f for f in os.listdir(tmpFldr) if os.path. \
		splitext(f)[1]=='.html']
		results = []
		hdr = False
		for rslts in lRslts:
			opndRstlsFl = open(os.path.join(tmpFldr,rslts)).read()
			soup = BeautifulSoup.RobustHTMLParser(opndRstlsFl)
			isHdr = True
			for tds in soup.findAll('tr'):
				if isHdr:
					hdr = '\t'.join([v.getText().encode('utf-8') for v \
					in tds.findAll('td')])
					isHdr = False
				else:
					vals = [str(v.getText()) for v in tds.findAll('td')]
					if vals[-1] != 'Non-Target':
						results.append('\t'.join(vals))
		results.sort()
		if hdr:
			results.insert(0,hdr)
		try:
			shutil.rmtree(tmpFldr)
		except:
			pass
		return results
	#
	def mthrdWrtTmp(qInJobs,qOutRslts,extnsnTrgtMnr,pthTmpFldr, \
		pthOutFldrTrgtMnr):
		"""
		Write a temporary file to run TargetMiner in a queue
		"""
		for dtNm in iter(qInJobs.get,'STOP'):			
			mrgRslTrgtMnrFl = ''.join([dtNm,extnsnTrgtMnr])
			mrgRslTrgtMnr = mrgRslts(dtNm,pthTmpFldr)
			mrgRslTrgtMnrFlNm = os.path.join(pthOutFldrTrgtMnr, \
			mrgRslTrgtMnrFl)
			mrgRslTrgtMnrFlNm = open(mrgRslTrgtMnrFlNm,'w')
			mrgRslTrgtMnrFlNm.write('\n'.join(mrgRslTrgtMnr))
			mrgRslTrgtMnrFlNm.close()
			qOutRslts.put(0)
	#
	assert pthTrgtMnr
	#Parse ran results
	sDtNms = set()
	for dtNm in aDtTrgtMnr:
		tmpFldr = os.path.join(pthTmpFldr,'%s.d'%dtNm)
		if os.path.exists(tmpFldr):
			sDtNms.add(dtNm)
	if sDtNms:
		cntTmpFl = len(sDtNms)
		qInJobs = Queue()
		while len(sDtNms)>0:
			dtNm = sDtNms.pop()
			qInJobs.put(dtNm)
		#write in multithread
		qOutRslts = Queue()
		for p in xrange(nThrds):
			Process(target=mthrdWrtTmp,args=(qInJobs,qOutRslts, \
			extnsnTrgtMnr,pthTmpFldr,pthOutFldrTrgtMnr)).start()
		for p in xrange(cntTmpFl):
			qOutRslts.get()
		for p in xrange(nThrds):
			qInJobs.put('STOP')		
	#Write those that haven't been written
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	sOutFlsTrgtMnr = set()
	for dtNm in aDtTrgtMnr:
		outFlNm = ''.join([dtNm,extnsnTrgtMnr])
		UTRFrmtdFlTrgtMnrNm = os.path.join(pthFldrUTRFrmtdFlsTrgtMnr, \
		''.join([dtNm,extnsnUTRTrgtMnr]))
		outFl = os.path.join(pthOutFldrTrgtMnr,outFlNm)
		sOutFlsTrgtMnr.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			wrtTmp(dtNm,UTRFrmtdFlTrgtMnrNm,cntMirnasGrps,
			mirnaFrmtdFlTrgtMnr,pthFldrMirnsFrmtdFlsTrgtMnr, \
			extnsnMirnaTrgtMnr,pthFldrPrwsFrmtdFlsTrgtMnr, \
			extnsnPrwsTrgtMnr,outFl,pthTmpFldr,tmpFl,pthTrgtMnr)
	#		
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:										
		#
		logFls = [os.path.join(pthTmpFldr,f) for f in \
		os.listdir(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsTrgtMnr


########################################################
#~ Run TargetScan in multithread
########################################################
def queueRunTrgtScn(aUTRFrmtdFls,mrnTrgtScnFrmtdFl,pthFldrUTRFrmtdFls, \
	pthFldrMirnsFrmtdFls,pthOutFldr,pthTmpFldr,pthTrgtScn,extnsnTrgtScn, \
	nPrlProcss):
	"""
	Input: aUTRFrmtdFls is an array of UTR TargetScan-formated file 
	names. mrnTrgtScnFrmtdFl is a miRNA TargetScan-formated file names. 
	pthFldrUTRFrmtdFls is a path to the folder with the UTR TargetScan-
	formated files. pthFldrMirnsFrmtdFls is a path to the folder with 
	the file in mrnTrgtScnFrmtdFl. pthOutFldr is a path to the folder to 
	save the TargetScan results. pthTmpFldr is a path to a temporal 
	folder to save the queue input files. pthTrgtScn is the path to the 
	TargetScan executable. extnsnTrgtScn is the extension of TargetScan 
	file results. nPrlProcss is the number of parallel processes to run 
	TargetScan. 
	Output: The results of targetScan written in the pthOutFldr. 
	sOutFlsTrgtScn is a the set of names for these files.
	"""
	#
	def wrtTmp(pthUTRFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthTrgtScn):
		"""
		Write a temporary file to run TargetScan in a queue
		"""
		fOpn = open(tmpFl,'w')
		fOpn.write('#!/bin/bash\n#\nexport PATH=%s:$PATH\n'%pthTrgtScn)
		fOpn.write('targetscan_70.pl %s %s %s\n'%(pthMirnsFrmtdFl, \
		pthUTRFrmtdFl,outFl))
		fOpn.close()
	#
	assert os.path.exists(pthTrgtScn)
	#
	sOutFlsTrgtScn = set()
	cntTmp = 0
	tmpFlNm	= NamedTemporaryFile().name[-9:]#with 6 random chars
	pthMirnsFrmtdFl = os.path.join(pthFldrMirnsFrmtdFls, \
	mrnTrgtScnFrmtdFl)
	for UTRFrmtdFl in aUTRFrmtdFls:
		UTRFrmtdFlPrfx = os.path.splitext(UTRFrmtdFl)[0]
		outFlNm = ''.join([UTRFrmtdFlPrfx,extnsnTrgtScn])
		outFl = os.path.join(pthOutFldr,outFlNm)
		sOutFlsTrgtScn.add(outFlNm)
		if not os.path.exists(outFl):
			cntTmp+=1
			tmpFl = os.path.join(pthTmpFldr,'%s_%s.sh'%(tmpFlNm,cntTmp))
			pthUTRFrmtdFl = os.path.join(pthFldrUTRFrmtdFls,UTRFrmtdFl)
			wrtTmp(pthUTRFrmtdFl,pthMirnsFrmtdFl,outFl,tmpFl,pthTrgtScn)
	#
	oufFlScrpt = os.path.join(pthTmpFldr,'%s_master.sh'%tmpFlNm)
	if cntTmp:
		oufFlScrpt = open(oufFlScrpt,'w')
		oufFlScrpt.write('#!/bin/bash -f\n\nsh %s_$SGE_TASK_ID.sh\n'% \
		tmpFlNm)
		oufFlScrpt.close()
		raise exceptions.CelleryExceptionObjct \
		('Go to path "%s"\nThen execute "qsub -t 1-%s -tc %s %s' \
		%(pthTmpFldr,cntTmp,nPrlProcss, \
		'-wd %s %s_master.sh"\n then run again this method'% \
		(pthTmpFldr,tmpFlNm)))
	else:
		logFls = [os.path.join(pthTmpFldr,f) for f in \
		os.listdir(pthTmpFldr) if f.find('.sh')>-1]
		for echfl in logFls:
			os.remove(echfl)
		return sOutFlsTrgtScn


########################################################
#~ Wrapper for running miRanda from fasta and MAF files and make 
#~ sql databases
########################################################
def runMiRanda(extnsnMAF,extnsnMiranda,extnsnUTRFsta,fstaDtFl, \
	fstasMirnaBsFl,nPrlProcss,nThrds,pthFldrFstaDtFl, \
	pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFstaOnlySppRef, \
	pthOutFldrMiranda,pthMiranda,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	sppRef,fOutDtSql,sqlFlMirnd,lAttrbts=['aMirnaMirndCnts', \
	'aMirnaMirndScrs','aMirnaMirndEngy']):
	"""
	Runs miRanda. 
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format MiRanda
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlMirnd):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlMirnd)
	else:
		#UTR formated files for MiRanda
		print 'Obtaining alignments...'
		alignments.proccsMAFnUTRsOnlySppRef(lDtObjcts,pthToMAF,extnsnMAF, \
		extnsnUTRFsta,pthFldrUTRFrmtdFlsFstaOnlySppRef,sppRef,pthToUCSCKEGG, \
		nThrds)
		aFstfrmtdFls = array(os.listdir(pthFldrUTRFrmtdFlsFstaOnlySppRef))
		#run miRanda
		print 'Running miRanda...'	
		sOutFlsRNAhybrd = queueRunRNAhybrd(aFstfrmtdFls,fstasMirnaBsFl, \
		pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthOutFldrMiranda, \
		pthTmpFldr,pthMiranda,extnsnMiranda,nPrlProcss)
		#parse RNAhybrid results	
		print 'Parsing miRanda results...'
		parsers.prsMirandaRslts(pthOutFldrMiranda,srtdMirnaFms, \
		extnsnMiranda,lDtObjcts)
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlMirnd)
	return lDtObjcts


########################################################
#~ Wrapper for running mirMap from fasta and MAF files and make 
#~ sql databases
########################################################
def runMirMap(extnsnMAF,extnsnMirMap,extnsnUTRFsta,fstaDtFl, \
	fstasMirnaBsFl,nPrlProcss,nThrds,phastModFl,phastTree, \
	pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFsta, \
	pthMirMapLib,pthOutFldrMirMap,pthPython,pthPythonLib,pthPhast, \
	pthTmpFldr,pthToMAF,pthToUCSCKEGG,sppRef,fOutDtSql,sqlFlMirMap, \
	lAttrbts=['aMirnaMirMapScrs','aMirnaMirMapCnts']):
	"""
	Runs miRmap. 
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format miRMap
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlMirMap):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlMirMap)
	else:
		#UTR formated files for making lObjects
		print 'Obtaining alignments...'
		alignments.proccsMAFnUTRsMirMap(lDtObjcts,pthToMAF,extnsnMAF, \
		extnsnUTRFsta,pthFldrUTRFrmtdFlsFsta,sppRef,pthToUCSCKEGG, \
		nThrds)
		aFstfrmtdFls = array(os.listdir(pthFldrUTRFrmtdFlsFsta))
		#run miRMap
		print 'Running mirMap...'
		sOutFlsMirMap = queueMirMap(aFstfrmtdFls,fstasMirnaBsFl, \
		pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthOutFldrMirMap,pthTmpFldr, \
		pthPython,pthPythonLib,pthPhast,pthMirMapLib,extnsnMirMap, \
		extnsnUTRFsta,phastModFl,phastTree,nPrlProcss)
		#parse miRMap results	
		print 'Parsing mirMap results...'
		parsers.prsMirMapRslts(pthOutFldrMirMap,srtdMirnaFms, \
		extnsnMirMap,lDtObjct)
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlMirMap)
	return lDtObjcts


########################################################
#~ Wrapper for running PITA from fasta and MAF files and make sql 
#~ databases
########################################################
def runPITA(extnsnMAF,extnsnUTRFsta,extnsnPITA,fstaDtFl,fstasMirnaBsFl, \
	nPrlProcss,nThrds,pthFldrFstaDtFl,pthFldrFstaMirnsFls, \
	pthFldrUTRFrmtdFlsFstaOnlySppRef,pthOutFldrPITA,pthPITA,pthTmpFldr, \
	pthToMAF,pthToUCSCKEGG,sppRef,fOutDtSql,sqlFlPITA, \
	lAttrbts=['aMirnaPITACnts','aMirnaPITAScrs']):
	"""
	Runs PITA. 
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format PITA
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	if os.path.exists(sqlFlPITA):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlPITA)
	else:
		#UTR formated files for making lObjects
		print 'Obtaining alignments...'
		alignments.proccsMAFnUTRsOnlySppRef(lDtObjcts,pthToMAF, \
		extnsnMAF,extnsnUTRFsta,pthFldrUTRFrmtdFlsFstaOnlySppRef,sppRef, \
		pthToUCSCKEGG,nThrds)
		aFstfrmtdFls = array(os.listdir \
		(pthFldrUTRFrmtdFlsFstaOnlySppRef))
		#run PITA
		print 'Running PITA...'	
		sOutFlsPITA = queueRunRNAhybrd(aFstfrmtdFls,fstasMirnaBsFl, \
		pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthOutFldrPITA, \
		pthTmpFldr,pthPITA,extnsnPITA,nPrlProcss)
		#import results from PITA and set pointers
		parsers.prsPITARslts(pthOutFldrPITA,srtdMirnaFms,extnsnPITA, \
		lDtObjct)
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlPITA)
	return lDtObjcts

########################################################
#~ Wrapper for running RNAhybrid from fasta and MAF files and make 
#~ sql databases
########################################################
def runRNAhybrid(extnsnMAF,extnsnRNAhybrd,extnsnUTRFsta,fstaDtFl, \
	fstasMirnaBsFl,nPrlProcss,nThrds,pthFldrFstaDtFl,pthFldrFstaMirnsFls, \
	pthFldrUTRFrmtdFlsFstaOnlySppRef,pthOutFldrRNAhybrd,pthRNAhybrd, \
	pthTmpFldr,pthToMAF,pthToUCSCKEGG,sppRef,fOutDtSql,sqlFlRNAhybrid, \
	lAttrbts=['aMirnaRNAhCnts','aMirnaRNAhEngy']):
	"""
	Runs RNAhybrid.
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format RNAhybrid
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlRNAhybrid):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlRNAhybrid)
	else:
		#UTR formated files for RNAhybrid
		print 'Obtaining alignments...'
		alignments.proccsMAFnUTRsOnlySppRef(lDtObjcts,pthToMAF,extnsnMAF, \
		extnsnUTRFsta,pthFldrUTRFrmtdFlsFstaOnlySppRef,sppRef, \
		pthToUCSCKEGG,nThrds)
		aFstfrmtdFls = array(os.listdir \
		(pthFldrUTRFrmtdFlsFstaOnlySppRef))
		#run RNAhybrid
		print 'Running RNAhybrid...'	
		sOutFlsRNAhybrd = queueRunRNAhybrd(aFstfrmtdFls,fstasMirnaBsFl, \
		pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthOutFldrRNAhybrd, \
		pthTmpFldr,pthRNAhybrd,extnsnRNAhybrd,nPrlProcss)
		#parse RNAhybrid results	
		print 'Parsing RNAhybrid results...'
		parsers.prsRNAhybrdRslts(pthOutFldrRNAhybrd,srtdMirnaFms, \
		extnsnRNAhybrd,lDtObjcts)
		#make sql file of MRE predictions
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlRNAhybrid)
	return lDtObjcts


########################################################
#~ Wrapper for running SVMicro from fasta and MAF files and make 
#~ sql databases
########################################################
def runSVMicro(extnsnMAF,extnsnPhstCns,extnsnSVMicro,fstaDtFl, \
	fstasMirnaBsFl,nPrlProcss,nThrds,pthFldrFstaDtFl, \
	pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsSVMicro,pthOutFldrSVMicro, \
	pthToPhstCns,pthSVMicro,pthTmpFldr,pthToMAF,pthToUCSCKEGG,sppRef, \
	fOutDtSql,sqlFlSVMicro,lAttrbts=['aMirnaSVMicroCnts', \
	'aMirnaSVMicroScrs']):
	"""
	Runs SVMicro. 
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format SVMicro
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlSVMicro):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlSVMicro)
	else:
		#UTR formated files for SVMicro
		print 'Obtaining alignments...'
		#format UTR sequences for SVMicro.
		alignments.procssUTRSVMicro(lDtObjcts,pthToPhstCns,extnsnPhstCns, \
		pthFldrUTRFrmtdFlsSVMicro,extnsnSVMicro,pthToMAF,extnsnMAF, \
		sppRef,pthToUCSCKEGG,nThrds)
		aUTRFrmtdFls = array(os.listdir(pthFldrUTRFrmtdFlsSVMicro))
		#run SVMicro
		print 'Running SVMicro...'	
		sOutFlsSVMicro = queueRunSVMicro(aUTRFrmtdFls,fstasMirnaBsFl, \
		pthFldrUTRFrmtdFlsSVMicro,pthFldrFstaMirnsFls,pthOutFldrSVMicro, \
		pthTmpFldr,pthSVMicro,extnsnSVMicro,nPrlProcss,nThrds)
		#import results from SVMicro and set pointers
		parsers.prsSVMicroRslts(pthOutFldrSVMicro,srtdMirnaFms, \
		extnsnSVMicro,lDtObjcts)
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlSVMicro)
	return lDtObjcts


########################################################
#~ Wrapper for running TargetMiner from fasta and MAF files and make 
#~ sql databases
########################################################
def runTargetMiner(extnsnMAF,extnsnMirnTrgtMnr,extnsnPrwsTrgtMnr, \
	extnsnTrgtMnr,extnsnUTRTrgtMnr,fstaDtFl,fstasMirnaBsFl, \
	mirnFrmtdFlTrgtMnr,nPrlProcss,nThrds,pthFldrFstaDtFl, \
	pthFldrFstaMirnsFls,pthFldrMirnsFrmtdFlsTrgtMnr, \
	pthFldrPrwsFrmtdFlsTrgtMnr,pthFldrUTRFrmtdFlsTrgtMnr, \
	pthOutFldrTrgtMnr,pthTrgtMnr,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	sppRef,fOutDtSql,sqlFlTrgtMnr,lAttrbts=['aMirnaTrgtMnrCnts']):
	"""
	Runs TargetMiner. 
	"""
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format TargetMiner
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,sqlFl=fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlTrgtMnr):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlTrgtMnr)
	else:
		#UTR formated files for making lObjects
		print 'Obtaining alignments...'
		alignments.procssUTRTrgtMnr(lDtObjcts,pthToMAF,extnsnMAF, \
		extnsnUTRTrgtMnr,pthFldrUTRFrmtdFlsTrgtMnr,sppRef,pthToUCSCKEGG, \
		nThrds)
		#format mirna file for TargetMiner
		lMirnasGrps = formats.procssMirnTrgtMnr(pthFldrFstaMirnsFls, \
		fstasMirnaBsFl,mirnFrmtdFlTrgtMnr,pthFldrMirnsFrmtdFlsTrgtMnr, \
		extnsnMirnTrgtMnr)
		cntMirnasGrps = len(lMirnasGrps)
		#make pairwise table of comparisons for TargetMiner	
		formats.mkPrwsFlTrgtMnr(lDtObjcts,lMirnasGrps, \
		pthFldrPrwsFrmtdFlsTrgtMnr,extnsnPrwsTrgtMnr,nThrds)
		aDtTrgtMnr = array([os.path.splitext(f)[0] for f in \
		os.listdir(pthFldrUTRFrmtdFlsTrgtMnr)])
		#run TargetMiner
		print 'Running TargetMiner...'	
		sOutFlsTrgtMnr = queueRunTrgtMnr(aDtTrgtMnr,cntMirnasGrps, \
		pthFldrUTRFrmtdFlsTrgtMnr,extnsnUTRTrgtMnr,mirnFrmtdFlTrgtMnr, \
		pthFldrMirnsFrmtdFlsTrgtMnr,extnsnMirnaTrgtMnr, \
		pthFldrPrwsFrmtdFlsTrgtMnr,extnsnPrwsTrgtMnr,pthOutFldrTrgtMnr, \
		pthTmpFldr,pthTrgtMnr,extnsnTrgtMnr,nPrlProcss,nThrds)
		#parse TargetMiner results
		prsTrgtMnrRslts(pthOutFldrTrgtMnr,srtdMirnaFms,extnsnTrgtMnr, \
		lDtObjct)
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlTrgtMnr)
	return lDtObjcts


########################################################
#~ Wrapper for running TargetScan from fasta and MAF files and make 
#~ sql databases
########################################################
def runTrgtScn(extnsnMAF,extnsnMirnaTrgtScn,extnsnTrgtScn, \
	extnsnUTRTrgtScn,fstaDtFl,fstasMirnaBsFl,mirbaseSpp, \
	nPrlProcss,nThrds,pthFldrFstaDtFl,pthFldrFstaMirnsFls, \
	pthFldrMirnsFrmtdFlsTrgtScn,pthFldrUTRFrmtdFlsTrgtScn, \
	pthOutFldrTrgtScn,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	pthTrgtScn,sppRef,sppMltzCd,fOutDtSql,sqlFlTrgtScn, \
	lAttrbts=['aMirnaCnts']):
	"""
	Runs TargetScan.
	"""	
	#sorted families of miRNAs including spp (i.e. miR-1178-5p)
	pthFlMirnaFms = os.path.join(pthFldrFstaMirnsFls,fstasMirnaBsFl)
	srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
	in open(pthFlMirnaFms).read().split('>') if l.strip()])
	#Execute for Genes
	print '\nExecuting analysis...'
	#format TargetScan
	if os.path.exists(fOutDtSql):
		sGnIntrst = set([s.split('|')[0][1:] for s in open(pthFlMirnaFms) \
		.readlines() if s[0]=='>'])
		lDtObjcts = __object__.mkArrySqlObjcts(fOutDtSql,srtdMirnaFms, \
		sGnIntrst)
	else:
		lDtObjcts = __object__.mkArryGnObjcts(fstaDtFl,srtdMirnaFms, \
		pthFldrFstaDtFl,fOutDtSql)
	#retrieve info for MREs
	if os.path.exists(sqlFlTrgtScn):
		#Retrieve info from SQL file
		for attrbt in lAttrbts:
			__object__.rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFlTrgtScn)
	else:
		#UTR formated files for TrgtScn. 
		print 'Obtaining alignments...'
		alignments.proccsMAFnUTRsTrgtScn(lDtObjcts,pthToMAF,extnsnMAF, \
		extnsnUTRTrgtScn,pthFldrUTRFrmtdFlsTrgtScn,sppRef,pthToUCSCKEGG, \
		nThrds)
		aUTRFrmtdFlsTrgtScnDt = array(os.listdir \
		(pthFldrUTRFrmtdFlsTrgtScn))
		#miRNA formated files for TrgtScn
		mrnTrgtScnFrmtdFl = ''.join([os.path.splitext(fstasMirnaBsFl)[0], \
		extnsnMirnaTrgtScn])
		tgrtScnFlmrnFrmtd = os.path.join(pthFldrMirnsFrmtdFlsTrgtScn, \
		mrnTrgtScnFrmtdFl)
		formats.mkArryFlsMirnaTrgtScn(pthFlMirnaFms,tgrtScnFlmrnFrmtd, \
		mirbaseSpp,sppMltzCd)
		#run TargetScan
		print 'Running TargetScan...'
		sOutFlsTrgtScn = queueRunTrgtScn(aUTRFrmtdFlsTrgtScnDt, \
		mrnTrgtScnFrmtdFl,pthFldrUTRFrmtdFlsTrgtScn, \
		pthFldrMirnsFrmtdFlsTrgtScn,pthOutFldrTrgtScn,pthTmpFldr, \
		pthTrgtScn,extnsnTrgtScn,nPrlProcss)
		print 'Parsing TargetScan results...'
		#import results from TargetScan and update objects...
		parsers.prsTrgtScnRslts(pthOutFldrTrgtScn,sppMltzCd,srtdMirnaFms, \
		extnsnTrgtScn,lDtObjcts)
		#make sql file of MRE predictions
		for attrbt in lAttrbts:
			__object__.mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFlTrgtScn)
	return lDtObjcts




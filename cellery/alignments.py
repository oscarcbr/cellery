#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  alignments.py part of cellery (ceRNAs linking inference)
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
Library with methods for retrieving alignments for genes/lncrnas or 
regions in the genome in different formats. For fasta format choose 
mirMap format.
"""

########################################################
#~ Import external libraries
########################################################
from multiprocessing import Queue, Process
from numpy import array, zeros, float32
from string import upper, lower
from Bio import Seq

import gzip
import os


########################################################
#~ Append sequences to a dictionary of spp names as keys and sequences 
# as values.
########################################################
def apndSeqs(dNmSeqs, cStrtEnd, sppRef, blckToAppnd, sUcscSpp):
	"""Input: dNmSeqs is a dictionary of stored {names:sequences}, 
	cStrtEnd tuple of the genome start and end of the current sequence 
	of reference, sppRef is the species of reference, blckToAppnd is 
	the block of information from the maf file. sUcscSpp is the set of 
	all species.
	Output: dNmSeqs is the updated dictionary of stored 
	{names:sequences} with the updated end, cStrtEnd is the updated 
	genome start and end of the stored sequences.
	"""
	cStrt, cEnd = cStrtEnd
	gapToAdd = False
	refSeq = False
	cntSpp = -1
	sSpp = set([ spp for spp in sUcscSpp ])#copy the set of spp
	for echl in blckToAppnd:
		if echl.strip() and echl[0] == 's':#sequence line. 
			cntSpp += 1
			typ, spp, gnmStrt, end, sense, othr, seq = echl.split()
			spp = spp.split('.')[0]
			if spp == sppRef:
				assert cntSpp == 0#test assume ref is 1st.
				gnmStrt, end = int(gnmStrt), int(end)
				end += gnmStrt
				cStrtEnd = (cStrt, end)
				if gnmStrt != cEnd:# in case there is missed sequences
					gapToAdd = 'N' * (gnmStrt - cEnd)
					seq = ''.join([gapToAdd, seq])
				refSeq = seq
			elif gapToAdd:#in case the sequences need extra gaps
				seq = ''.join([gapToAdd, seq])
			cSeq = dNmSeqs[spp]
			seq = ''.join([cSeq, seq])
			dNmSeqs[spp] = upper(seq)
			sSpp.remove(spp)
	#add sequences for the missing species
	if sSpp and refSeq:
		fllGapSeq = '-' * len(refSeq)
		while sSpp:
			spp = sSpp.pop()
			cSeq = dNmSeqs[spp]
			seq = ''.join([cSeq, fllGapSeq])
			dNmSeqs[spp] = seq
	return dNmSeqs, cStrtEnd


########################################################
#~ Merge sequences in a dictionary of spp names as keys and sequences 
# as values.
########################################################
def assmblGnExns(dIntrvldNmSeqs, sppRef, strnd, gn, TforU=False):
	"""Input: dIntrvldNmSeqs is a dictionary of intervals {spp:seq}.
	sppRef is the species of reference, strnd is the strnd of the gene, 
	if TforU is True all Ts are replaced by Us.
	Output: dMrgdNmSeqs is dictionary of {spp:seq} with the sequences
	merged. lErrs includes the list of gene names and intervals that
	only retrieved gaps for the reference.
	"""
	lErrs = []
	#make a set of spp in at least 1 alignment
	sSpp = set()
	for intrvl in dIntrvldNmSeqs.keys():
		sSpp.update(set(dIntrvldNmSeqs[intrvl].keys()))
	#append alignments
	lIntrvls = dIntrvldNmSeqs.keys()
	lIntrvls.sort()
	dMrgdNmSeqs = dict([ (spp, []) for spp in sSpp ])
	while lIntrvls:
		intrvl = lIntrvls.pop(0)
		dNmSeqs = dIntrvldNmSeqs.pop(intrvl)
		seq = dNmSeqs[sppRef]
		if set(seq) == {'-'} or set(seq) == {'N'}:#not only gaps
			lErrs.append('\t'.join([gn, str(intrvl[0]), \
			str(intrvl[1])]))
		emptyAlgnmnt = '-' * len(seq)
		for spp in sSpp:
			if dNmSeqs.has_key(spp):
				seq = dNmSeqs[spp]
				if TforU:
					seq = seq.replace('T', 'U')
				dMrgdNmSeqs[spp].append(seq)
			else:
				dMrgdNmSeqs[spp].append(emptyAlgnmnt)
	#merge alignments
	if strnd == '1':
		dMrgdNmSeqs = dict([ (spp, ''.join(dMrgdNmSeqs[spp])) for spp \
		in sSpp ])
	else:
		dMrgdNmSeqs = dict([ (spp, Seq.reverse_complement(''.join \
		(dMrgdNmSeqs[spp]))) for spp in sSpp ])
	return dMrgdNmSeqs, lErrs


########################################################
#~ Make a dictionary of chromosomes as keys and lists of start, end, 
# and names in genes as values
########################################################
def mkdChrlStrtEndGnNmStrnd(lDtObjcts, pthFldrUTRFrmtdFls, extnsnUTR):
	"""Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "chr", and "strnd" for the genes or 
	lncrnas. pthFldrUTRFrmtdFls is the path to the folder to output the 
	gene/lncrna sequences in different formats for each chromosome. 
	extnsnUTR is the extension of the output UTR files. 
	Output: dChrlStrtEndGnNm is a dictionary of chromosomes as keys and 
	as values lists of start, end, and name. dGnNmStrnd is a dictionary 
	of {gene name:strand}.
	"""
	dChrlStrtEndGnNm = {}
	dGnNmStrnd = {}
	for dt in lDtObjcts:
		name = dt.name
		outFlFl = os.path.join(pthFldrUTRFrmtdFls, ''.join([name, \
		extnsnUTR]))
		if not os.path.exists(outFlFl) or not os.path.getsize(outFlFl):
			#skip those in path
			chr = dt.chr
			strnd = dt.strnd
			dGnNmStrnd[name] = strnd
			for intrvl in dt.aIntrvls:
				strt, end = intrvl
				if dChrlStrtEndGnNm.has_key(chr):
					dChrlStrtEndGnNm[chr].append((strt, end, name))
				else:
					dChrlStrtEndGnNm[chr] = [(strt, end, name)]
	return dChrlStrtEndGnNm, dGnNmStrnd


########################################################
#~ Multithread processing of retrieving sequence alignments for each
# chromosome.
########################################################
def mthrdRtrvMAFs(qInJobs, qOutRslts, sUcscSpp, sppRef, dGnNmStrnd, \
	pthFldrUTRFrmtdFls, extnsnUTR, wrtAlgnmnt, TforU=False, \
	dUcscTOMultizSpp=False):
	"""Input: qInJobs is a queue with (fMAF,lStrtEndGnNm). fMAF is the 
	full path to the MAF file. lStrtEndGnNm is a list of start, end, 
	and gene names. sUcscSpp is a set of species to be included in the 
	alignment (must include all in the MAF file). sppRef is a reference 
	species. dGnNmStrnd is a dictionary of {gene name:strand}. 
	pthFldrUTRFrmtdFls is the outfolder to save the alignments, 
	extnsnUTR is the extension of alignments output files. qOutRslts is 
	a queue to output the list of errors. Optionally, dUcscTOMultizSpp 
	is dictionary of multiz species and TargetScan codes.
	Output: lErrs is a list for which part of the alignment for 
	the reference gene was found to be empty. lExcd is a list of 
	excluded intervals outside MAF files.
	"""
	for fMAF, lStrtEndGnNm in iter(qInJobs.get, 'STOP'):
		lErrs, lExcd = runCore(fMAF, lStrtEndGnNm, sUcscSpp, sppRef, \
		dGnNmStrnd, pthFldrUTRFrmtdFls, extnsnUTR, wrtAlgnmnt, TforU, \
		dUcscTOMultizSpp=dUcscTOMultizSpp)
		qOutRslts.put((lErrs, lExcd))


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest genome-wide in fasta/mirMap format.
########################################################
def proccsMAFnUTRsMirMap(lDtObjcts, pthToMAF, extnsnMAF, \
	extnsnUTRMirMap, pthFldrUTRFrmtdFlsMirMap, sppRef, pthToUCSCKEGG, \
	nThrds=5):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "chr", and "strnd" for the genes 
	or lncrnas. pthToMAF is a path to the folder with the input MAF 
	files for each chromosome. extnsnMAF is the extension for the MAF 
	files in the pthToMAF path. pthFldrUTRFrmtdFlsMirMap is the path to 
	the folder to output the gene/lncrna sequences in miRMap/fasta 
	format for each chromosome. sppRef is the UCSC name of the species 
	(i.e. hg19). nThrds is the number of parallel processes to run 
	mirMap. extnsnUTRMirMap is the extension of the output UTR files. 
	pthToUCSCKEGG is the path to the files relating all UCSC codes to
	other codes.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", and "seq" in the same order as the sorted names of 
	for the genes/lncrnas of interest in lDtObjcts. It also writes the
	miRMap/fasta  formated sequence files in pthFldrUTRFrmtdFlsMirMap. 
	Only those genes/lncrnas with information are going to be included 
	in the output. A list of errors in printed (if any).
	NOTE: A high number of threads requires a higher memory usage.
	"""
	#
	def wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFlsMirMap, \
		extnsnUTRMirMap, sppRef):
		""" 
		Input: dGendMrgdNmSeqs is dictionary of {gene:{spp:seq}} 
		with the sequences merged, pthFldrUTRFrmtdFlsMirMap is the 
		output folder to save the alignments, extnsnUTRMirMap is
		the extension of the alignment file, sppRef is the reference 
		species. 
		This method will write the alignment sequences into the output
		folder in miRMap/fasta format.
		"""
		sGnNms = set(dGendMrgdNmSeqs.keys())#set of gene names
		while sGnNms:
			name = sGnNms.pop()
			outFlMirMapFl = os.path.join(pthFldrUTRFrmtdFlsMirMap, \
			''.join([name, extnsnUTRMirMap]))
			dSppAlgndSeqs = dGendMrgdNmSeqs.pop(name)
			refSeq = dSppAlgndSeqs.pop(sppRef)
			lFrtmdUTR = ['>%s\n%s' % (sppRef, refSeq)]
			sSpp = set(dSppAlgndSeqs.keys())#set of species
			while sSpp:
				spp = sSpp.pop()
				seq = dSppAlgndSeqs.pop(spp)
				lFrtmdUTR.append('>%s\n%s' % (spp, seq))
			#write file
			outFlMirMapFl = open(outFlMirMapFl, 'w')
			outFlMirMapFl.write('\n'.join(lFrtmdUTR))
			outFlMirMapFl.close()
		return 0
	#--------------------------
	#Start variables to test
	sUcscSpp = set([ l.split()[0] for l in open(pthToUCSCKEGG).read(). \
	splitlines() if l.strip() ])
	#Retrieve dictionary of sorted information by chromosome
	dChrlStrtEndGnNm, dGnNmStrnd = mkdChrlStrtEndGnNmStrnd(lDtObjcts, \
	pthFldrUTRFrmtdFlsMirMap, extnsnUTRMirMap)
	#Clear the info for chromosomes not present in pthToMAF
	sChrs = set(dChrlStrtEndGnNm.keys())
	sChrsinMAFpth = set([ f.split(extnsnMAF)[0][3:] for f in os.listdir \
	(pthToMAF) if f.find(extnsnMAF) > -1 ])
	sChrAbsnt = sChrs.difference(sChrsinMAFpth)
	while sChrAbsnt:
		eChr = sChrAbsnt.pop()
		inf = dChrlStrtEndGnNm.pop(eChr)
		del inf
	#Retrieve alignments
	sChrsinMAFpth = sChrsinMAFpth.intersection(sChrs)
	lenChrs = len(sChrsinMAFpth)
	lErrsGlbl, lExcdGlbl = [], []
	qInJobs = Queue()
	qOutRslts = Queue()
	cnt = 0
	while sChrsinMAFpth:
		cnt += 1
		eChr = sChrsinMAFpth.pop()
		print ''.join(['\t...retrieving alignments for fasta/miRMap from', 
		' UTRs in chromosome %s, number %s out of %s' % (eChr, cnt, \
		lenChrs)])
		fMAF = os.path.join(pthToMAF, 'chr%s.maf.gz' % eChr)
		lStrtEndGnNm = dChrlStrtEndGnNm.pop(eChr)
		qInJobs.put((fMAF, lStrtEndGnNm))
	for thrd in xrange(nThrds):
		Process(target=mthrdRtrvMAFs, args=(qInJobs,
		 qOutRslts,
		 sUcscSpp,
		 sppRef,
		 dGnNmStrnd,
		 pthFldrUTRFrmtdFlsMirMap,
		 extnsnUTRMirMap,
		 wrtAlgnmnt)).start()
	cnt = 0
	for eChr in xrange(lenChrs):
		cnt += 1
		lErrs, lExcd = qOutRslts.get()
		lErrsGlbl.extend(lErrs)
		lExcdGlbl.extend(lExcd)
		print '\tRetrieving results for chromosome number %s' % cnt
	for thrd in xrange(nThrds):
		qInJobs.put('STOP')
	#print errors
	if lErrsGlbl:
		print '\tThe following sequence(s) for gene(s) and position(s) were emptied:'
		print '\t%s' % '\n\t'.join(lErrsGlbl)
		print '\tAlignments with gaps in these positions were retrieved.'
	if lExcdGlbl:
		print '\tThe following interval(s) for gene(s) and position(s) were out of MAF:'
		print '\t%s' % '\n\t'.join(lExcdGlbl)
		print '\tAlignments for the interval(s) were not retrieved.'
	return 0


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest genome-wide in fasta/RNAhybrid format.
########################################################
def proccsMAFnUTRsOnlySppRef(lDtObjcts, pthToMAF, extnsnMAF, \
	extnsnUTROnlySppRef, pthFldrUTRFrmtdFlsOnlySppRef, sppRef, \
	pthToUCSCKEGG, nThrds=5):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "chr", and "strnd" for the genes 
	or lncrnas. pthToMAF is a path to the folder with the input MAF 
	files for each chromosome. extnsnMAF is the extension for the MAF 
	files in the pthToMAF path. pthFldrUTRFrmtdFlsOnlySppRef is the path 
	to the folder to output the gene/lncrna sequences in OnlySppRef/
	fasta format for each chromosome. sppRef is the UCSC name of the 
	species (i.e. hg19). nThrds is the number of parallel processes to 
	run OnlySppRef. extnsnUTROnlySppRef is the extension of the output 
	UTR files. pthToUCSCKEGG is the path to the files relating all UCSC 
	codes to other codes. 
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", and "seq" in the same order as the sorted names of 
	for the genes/lncrnas of interest in lDtObjcts. It also writes the
	OnlySppRef/fasta  formated sequence files in 
	pthFldrUTRFrmtdFlsOnlySppRef. Only those genes/lncrnas with 
	information are going to be included in the output. A list of errors 
	in printed (if any).
	NOTE: A high number of threads requires a higher memory usage.
	NOTE: Only the sequence	for the reference species is going to be 
	returned.
	"""
	#
	def wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFlsOnlySppRef, \
		extnsnUTROnlySppRef, sppRef):
		""" 
		Input: dGendMrgdNmSeqs is dictionary of {gene:{spp:seq}} 
		with the sequences merged, pthFldrUTRFrmtdFlsOnlySppRef is the 
		output folder to save the alignments, extnsnUTROnlySppRef is
		the extension of the alignment file, sppRef is the reference 
		species. 
		This method will write the alignment sequences into the output
		folder in OnlySppRef/fasta format.
		"""
		sGnNms = set(dGendMrgdNmSeqs.keys())#set of gene names
		while sGnNms:
			name = sGnNms.pop()
			outFlOnlySppRefFl = os.path.join(pthFldrUTRFrmtdFlsOnlySppRef, \
			''.join([name, extnsnUTROnlySppRef]))
			dSppAlgndSeqs = dGendMrgdNmSeqs.pop(name)
			refSeq = dSppAlgndSeqs.pop(sppRef).replace('-','')
			lFrtmdUTR = ['>%s\n%s' % (name, refSeq)]
			#write file
			outFlOnlySppRefFl = open(outFlOnlySppRefFl, 'w')
			outFlOnlySppRefFl.write('\n'.join(lFrtmdUTR))
			outFlOnlySppRefFl.close()
		return 0
	#--------------------------
	#Start variables to test
	sUcscSpp = set([ l.split()[0] for l in open(pthToUCSCKEGG).read(). \
	splitlines() if l.strip() ])
	#Retrieve dictionary of sorted information by chromosome
	dChrlStrtEndGnNm, dGnNmStrnd = mkdChrlStrtEndGnNmStrnd(lDtObjcts, \
	pthFldrUTRFrmtdFlsOnlySppRef, extnsnUTROnlySppRef)
	#Clear the info for chromosomes not present in pthToMAF
	sChrs = set(dChrlStrtEndGnNm.keys())
	sChrsinMAFpth = set([ f.split(extnsnMAF)[0][3:] for f in os.listdir \
	(pthToMAF) if f.find(extnsnMAF) > -1 ])
	sChrAbsnt = sChrs.difference(sChrsinMAFpth)
	while sChrAbsnt:
		eChr = sChrAbsnt.pop()
		inf = dChrlStrtEndGnNm.pop(eChr)
		del inf
	#Retrieve alignments
	sChrsinMAFpth = sChrsinMAFpth.intersection(sChrs)
	lenChrs = len(sChrsinMAFpth)
	lErrsGlbl, lExcdGlbl = [], []
	qInJobs = Queue()
	qOutRslts = Queue()
	cnt = 0
	while sChrsinMAFpth:
		cnt += 1
		eChr = sChrsinMAFpth.pop()
		print ''.join(['\t...retrieving alignments for fasta/OnlySppRef', \
		' from UTRs in chromosome %s, number %s out of %s' % (eChr, cnt, \
		lenChrs)])
		fMAF = os.path.join(pthToMAF, 'chr%s.maf.gz' % eChr)
		lStrtEndGnNm = dChrlStrtEndGnNm.pop(eChr)
		qInJobs.put((fMAF, lStrtEndGnNm))
	for thrd in xrange(nThrds):
		Process(target=mthrdRtrvMAFs, args=(qInJobs,
		 qOutRslts,
		 sUcscSpp,
		 sppRef,
		 dGnNmStrnd,
		 pthFldrUTRFrmtdFlsOnlySppRef,
		 extnsnUTROnlySppRef,
		 wrtAlgnmnt)).start()
	cnt = 0
	for eChr in xrange(lenChrs):
		cnt += 1
		lErrs, lExcd = qOutRslts.get()
		lErrsGlbl.extend(lErrs)
		lExcdGlbl.extend(lExcd)
		print '\tRetrieving results for chromosome number %s' % cnt
	for thrd in xrange(nThrds):
		qInJobs.put('STOP')
	#print errors
	if lErrsGlbl:
		print '\tThe following sequence(s) for gene(s) and position(s) were emptied:'
		print '\t%s' % '\n\t'.join(lErrsGlbl)
		print '\tAlignments with gaps in these positions were retrieved.'
	if lExcdGlbl:
		print '\tThe following interval(s) for gene(s) and position(s) were out of MAF:'
		print '\t%s' % '\n\t'.join(lExcdGlbl)
		print '\tAlignments for the interval(s) were not retrieved.'
	return 0


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest genome-wide in TargetScan format.
########################################################
def proccsMAFnUTRsTrgtScn(lDtObjcts,pthToMAF,extnsnMAF,extnsnUTRTrgtScn, \
	pthFldrUTRFrmtdFlsTrgtScn, sppRef, pthToUCSCKEGG, nThrds=5):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "chr", and "strnd" for the genes 
	or lncrnas. pthToMAF is the path to the folder with the input MAF 
	files for each chromosome. extnsnMAF is the extension for the MAF 
	files in the pthToMAF path. pthFldrUTRFrmtdFlsTrgtScn is the path to 
	the folder to output the gene/lncrna sequences in  fomra for each 
	chromosome. sppRef is the UCSC name of the species (i.e. hg19).  
	nThrds is the number of parallel processes to run TargetScan. 
	extnsnUTRTrgtScn is the extension of the output UTR files. 
	pthToUCSCKEGG is the path to the files relating all UCSC codes to
	other codes.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", and "seq" in the same order as the sorted names of 
	for the genes/lncrnas of interest in lDtObjcts. It also writes the
	TargetScan formated sequence files in pthFldrUTRFrmtdFlsTrgtScn. 
	Only those genes/lncrnas with information are going to be included 
	in the output. A list of errors in printed (if any).
	NOTE: A high number of threads requires a higher memory usage.
	"""
	#
	def wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFlsTrgtScn, \
		extnsnUTRTrgtScn, sppRef, dUcscTOMultizSpp):
		""" Input: dGendMrgdNmSeqs is dictionary of {gene:{spp:seq}} 
		with the sequences merged, pthFldrUTRFrmtdFlsTrgtScn is the 
		output folder to save the alignments, extnsnUTRTrgtScn is
		the extension of the alignment file, sppRef is the reference 
		species, dUcscTOMultizSpp is a dictionary of spp names and 
		multiz codes.
		This method will write the alignment sequences into the output
		folder in TargetScan format.
		"""
		sGnNms = set(dGendMrgdNmSeqs.keys())#set of gene names
		while sGnNms:
			name = sGnNms.pop()
			outFlTrgtScnFl = os.path.join(pthFldrUTRFrmtdFlsTrgtScn, \
			''.join([name, extnsnUTRTrgtScn]))
			dSppAlgndSeqs = dGendMrgdNmSeqs.pop(name)
			refSeq = dSppAlgndSeqs.pop(sppRef)
			lFrtmdUTR = ['\t'.join([name, dUcscTOMultizSpp[sppRef], \
			refSeq])]
			sSpp = set(dSppAlgndSeqs.keys())#set of species
			while sSpp:
				spp = sSpp.pop()
				seq = dSppAlgndSeqs.pop(spp)
				multizspp = dUcscTOMultizSpp[spp]
				lFrtmdUTR.append('\t'.join([name, multizspp, seq]))
			#write file
			outFlTrgtScnFl = open(outFlTrgtScnFl, 'w')
			outFlTrgtScnFl.write('\n'.join(lFrtmdUTR))
			outFlTrgtScnFl.close()
		return 0
	#--------------------------
	#Start variables to test
	dUcscTOMultizSpp = dict([ (l.split()[0], l.split()[-1]) for l in \
	open(pthToUCSCKEGG).read().splitlines() if l.strip() ])
	sUcscSpp = set(dUcscTOMultizSpp.keys())
	dChrlStrtEndGnNm, dGnNmStrnd = mkdChrlStrtEndGnNmStrnd(lDtObjcts, \
	pthFldrUTRFrmtdFlsTrgtScn, extnsnUTRTrgtScn)
	sChrs = set(dChrlStrtEndGnNm.keys())
	sChrsinMAFpth = set([ f.split(extnsnMAF)[0][3:] for f in \
	os.listdir(pthToMAF) if f.find(extnsnMAF) > -1 ])
	sChrAbsnt = sChrs.difference(sChrsinMAFpth)
	while sChrAbsnt:
		eChr = sChrAbsnt.pop()
		inf = dChrlStrtEndGnNm.pop(eChr)
		del inf
	#Retrieve alignments
	sChrsinMAFpth = sChrsinMAFpth.intersection(sChrs)
	lenChrs = len(sChrsinMAFpth)
	lErrsGlbl, lExcdGlbl = [], []
	qInJobs = Queue()
	qOutRslts = Queue()
	cnt = 0
	while sChrsinMAFpth:
		cnt += 1
		eChr = sChrsinMAFpth.pop()
		print ''.join(['\t...retrieving alignments for TargetScan from', \
		' UTRs in chromosome %s, number %s out of %s' % (eChr, cnt, \
		lenChrs)])
		fMAF = os.path.join(pthToMAF, 'chr%s.maf.gz' % eChr)
		lStrtEndGnNm = dChrlStrtEndGnNm.pop(eChr)
		qInJobs.put((fMAF, lStrtEndGnNm))
	TforU = True
	for thrd in xrange(nThrds):
		Process(target=mthrdRtrvMAFs, args=(qInJobs,
		 qOutRslts,
		 sUcscSpp,
		 sppRef,
		 dGnNmStrnd,
		 pthFldrUTRFrmtdFlsTrgtScn,
		 extnsnUTRTrgtScn,
		 wrtAlgnmnt,
		 TforU,
		 dUcscTOMultizSpp)).start()
	cnt = 0
	for eChr in xrange(lenChrs):
		cnt += 1
		lErrs, lExcd = qOutRslts.get()
		lErrsGlbl.extend(lErrs)
		lExcdGlbl.extend(lExcd)
		print '\tRetrieving results for chromosome number %s' % cnt
	for thrd in xrange(nThrds):
		qInJobs.put('STOP')
	#print errors
	if lErrsGlbl:
		print '\tThe following sequence(s) for gene(s) and position(s) were emptied:'
		print '\t%s' % '\n\t'.join(lErrsGlbl)
		print '\tAlignments with gaps in these positions were retrieved.'
	if lExcdGlbl:
		print '\tThe following interval(s) for gene(s) and position(s) were out of MAF:'
		print '\t%s' % '\n\t'.join(lExcdGlbl)
		print '\tAlignments for the interval(s) were not retrieved.'
	return 0


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest genome-wide in SVMicro format.
########################################################
def procssUTRSVMicro(lDtObjcts, pthToPhstCns, extnsnPhstCns, \
	pthFldrUTRFrmtdFlsSVMicro, extnsnSVMicro, pthToMAF, extnsnMAF, \
	sppRef, pthToUCSCKEGG, nThrds=5):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "eChr", and "strnd" for the genes 
	or lncrnas. pthToPhstCns is the path to the folder with the input 
	phastCons score files for each chromosome. extnsnPhstCns is the 
	extension for the phastCons score files in the pthToPhstCns path. 
	pthFldrUTRFrmtdFlsSVMicro is the path to the folder to output the 
	gene/lncrna sequences in SVMicro format. extnsnSVMicro is the 
	extension of the output UTR files (in SVMicro format). pthToMAF is 
	a path to the folder with the input MAF files for each chromosome. 
	extnsnMAF is the extension for the MAF files in the pthToMAF path.
	sppRef is the UCSC name of the species (i.e. hg19). pthToUCSCKEGG is 
	the path to the files relating all UCSC codes to other codes. nThrds 
	is the number of parallel processes to run TargetScan.
	Output:  SVMicro formatted UTR sequences in the 
	pthFldrUTRFrmtdFlsSVMicro folder. aUTRFrmtdFlsSVMicro is an array 
	with the sorted set of files for the SVMicro results.
	"""

	def runCoreSVMicro(lDtObjcts, fMAF, lStrtEndGnNm, sUcscSpp, sppRef, \
		dGnNmStrnd, eChr, pthToPhstCns, extnsnPhstCns, \
		pthFldrUTRFrmtdFlsSVMicro, extnsnSVMicro):
		""" 
		Input: lDtObjcts is a list with gene/lncrna objects with pointers  
		"pos", "name", "aMrdgIntrvls", "eChr", and "strnd" for the genes 
		or lncrnas. fMAF is the input MAF files for the chromosome of 
		interest. lStrtEndGnNm is a list with a sorted tuple of (start, 
		end,gene_name) for each gene in the chromosome of interest. 
		sUcscSpp is a set of species to be included in the alignment 
		(must include all in the MAF file). sppRef is a reference 
		species. dGnNmStrnd is a dictionary of {gene name:strand}. eChr
		is the name of the chromosome. pthToPhstCns is the path to 
		phastCons files. extnsnPhstCns is the extension of the phastCons 
		files. pthFldrUTRFrmtdFlsSVMicro is the outfolder to save the 
		alignments in SVMicro format, extnsnSVMicro is the extension of 
		alignments SVMicro output files.
		Output: sFlSVMicroNms is the set of SVMicro files created, lErrs
		is a list of errors produced when retrieving the sequences. 
		lExcd is a list of excluded intervals outside MAF files.
		"""
		srtdNtPos = set()
		for strt, end, gnNm in lStrtEndGnNm:
			srtdNtPos.update(set(range(strt, end)))
		#Return phastCons score
		srtdNtPos = sorted(srtdNtPos)
		srtdNtPos.reverse()
		flPhtConsNm = ''.join(['chr', eChr, extnsnPhstCns])
		flPhtCons = os.path.join(pthToPhstCns, flPhtConsNm)
		fPhtCons = gzip.open(flPhtCons, 'rb')
		cNtPos = srtdNtPos.pop()
		phtCons = 0
		dPosPhstCons = {}
		breaker = False#nested loop breaker
		for echl in fPhtCons:
			if echl[0] == 'f':
				pfx, pfx, curPos, pfx = echl.split()
				curPos = int(curPos.split('=')[1]) - 1
			else:
				if curPos > cNtPos:
					while curPos > cNtPos:
						dPosPhstCons[cNtPos] = float32(0)
						if len(srtdNtPos) > 0:
							cNtPos = srtdNtPos.pop()
						else:
							breaker = True
							break
				if curPos == cNtPos:
					phtCons = float32(echl.strip())
					dPosPhstCons[curPos] = phtCons
					if len(srtdNtPos) > 0:
						cNtPos = srtdNtPos.pop()
					else:
						breaker = True
				if breaker:
					break
				else:
					curPos += 1
		#Return sequences for the chromosome
		dGendMrgdNmSeqs, lErrs, lExcd = runCore(fMAF, lStrtEndGnNm, \
		sUcscSpp, sppRef, dGnNmStrnd, 0, 0, 0, TforU=False, \
		rtrndGendMrgdNmSeqs=True)
		sGnsInChr = set(dGendMrgdNmSeqs.keys())
		#format SVMicro
		for dt in lDtObjcts:
			name = dt.name
			if name in sGnsInChr:
				lFrmtdSVMicro = ['gene_id gene_sym\tmrna_acc\tsequence\tconsv_score']
				aIntrvls = dt.aIntrvls
				seq = dGendMrgdNmSeqs[name][sppRef].replace('-','')
				lPhstCons = []
				for strt, end in aIntrvls:
					lPhstCons.extend([ str(dPosPhstCons[nt]) for nt in \
					xrange(strt, end) ])
				if dt.strnd[0] == '-':
					lPhstCons.reverse()
				try:
					assert len(seq)==len(lPhstCons)
				except:
					print len(seq),len(lPhstCons),name
				lFrmtdSVMicro.append('\t'.join(['1',name,name,seq, \
				', '.join(lPhstCons)]))
				flSVMicroNm = ''.join([name, extnsnSVMicro])
				flSVMicro = os.path.join(pthFldrUTRFrmtdFlsSVMicro, \
				flSVMicroNm)
				flSVMicro = open(flSVMicro, 'w')
				flSVMicro.write('\n'.join(lFrmtdSVMicro))
				flSVMicro.close()
		return lErrs, lExcd
	#
	def mthrdFrmtSVMicro(qInJobs, qOutRslts, lDtObjcts, sUcscSpp, \
		sppRef, dGnNmStrnd, pthToPhstCns, extnsnPhstCns, \
		pthFldrUTRFrmtdFlsSVMicro, extnsnSVMicro):
		"""
		Input: qInJobs is an is a queue with (eChr,fMAF,lStrtEndGnNm) 
		values for the chromosome of interest, qOutRslts. lDtObjcts is a 
		list with gene/lncrna objects with pointers "pos", "name", 
		"aMrdgIntrvls", "eChr", and "strnd" for the genes or lncrnas. 
		sppRef is a reference species.  fMAF is the input MAF files for 
		the chromosome of interest. UcscSpp is a set of species to be 
		included in the alignment (must include all in the MAF file). 
		lStrtEndGnNm is a list with a sorted tuple of (start,end,
		gene_name) for each gene in the chromosome of interest. eChr is 
		the name of the chromosome. dGnNmStrnd is a dictionary of {gene 
		name:strand}. pthToPhstCns is the path to phastCons files. 
		extnsnPhstCns is the extension of the phastCons files. 
		pthFldrUTRFrmtdFlsSVMicro is the outfolder to save the 
		alignments in SVMicro format, extnsnSVMicro is the extension of 
		alignments SVMicro output files.
		Output: qOutRslts is the output queue with values for and lErrs
		added. sFlSVMicroNms is the set of SVMicro files created, lErrs
		is a list of errors produced when retrieving the sequences. 
		lExcd is a list of excluded intervals outside MAF files.
		"""
		for eChr, fMAF, lStrtEndGnNm in iter(qInJobs.get, 'STOP'):
			lErrs, lExcd = runCoreSVMicro(lDtObjcts, fMAF, lStrtEndGnNm, \
			sUcscSpp, sppRef, dGnNmStrnd, eChr, pthToPhstCns, \
			extnsnPhstCns, pthFldrUTRFrmtdFlsSVMicro, extnsnSVMicro)
			qOutRslts.put((lErrs, lExcd))
	#--------------------------
	#Start variables to retrieve alignments
	dUcscTOMultizSpp = dict([ (l.split()[0], l.split()[-1]) for l in \
	open(pthToUCSCKEGG).read().splitlines() if l.strip() ])
	sUcscSpp = set(dUcscTOMultizSpp.keys())
	#Retrieve dictionary of sorted information by chromosome
	dChrlStrtEndGnNm, dGnNmStrnd = mkdChrlStrtEndGnNmStrnd(lDtObjcts, \
	pthFldrUTRFrmtdFlsSVMicro, extnsnSVMicro)
	#Clear the info for chromosomes not present in pthToMAF
	sChrs = set(dChrlStrtEndGnNm.keys())#set of total chromosomes
	sChrsinMAFpth = set([ f.split(extnsnMAF)[0][3:] for f in \
	os.listdir(pthToMAF) if f.find(extnsnMAF) > -1 ])
	sChrAbsnt = sChrs.difference(sChrsinMAFpth)#real
	while sChrAbsnt:
		eChr = sChrAbsnt.pop()
		inf = dChrlStrtEndGnNm.pop(eChr)
		del inf
	#Retrieve alignments
	sChrsinMAFpth = sChrsinMAFpth.intersection(sChrs)
	lenChrs = len(sChrsinMAFpth)
	lErrsGlbl, lExcdGlbl = [], []
	qInJobs = Queue()
	qOutRslts = Queue()
	cnt = 0
	while sChrsinMAFpth:
		cnt += 1
		eChr = sChrsinMAFpth.pop()
		print ''.join(['\t...retrieving alignments for SVMicro from', 
		' UTRs in chromosome %s, number %s out of %s' % (eChr, cnt, \
		lenChrs)])
		fMAF = os.path.join(pthToMAF, 'chr%s.maf.gz' % eChr)
		lStrtEndGnNm = dChrlStrtEndGnNm.pop(eChr)
		qInJobs.put((eChr, fMAF, lStrtEndGnNm))
	for thrd in xrange(nThrds):
		Process(target=mthrdFrmtSVMicro, args=(qInJobs,qOutRslts, \
		 lDtObjcts,sUcscSpp,sppRef,dGnNmStrnd,pthToPhstCns, \
		 extnsnPhstCns,pthFldrUTRFrmtdFlsSVMicro,extnsnSVMicro)).start()
	cnt = 0
	for eChr in xrange(lenChrs):
		cnt += 1
		lErrs, lExcd = qOutRslts.get()
		lErrsGlbl.extend(lErrs)
		lExcdGlbl.extend(lExcd)
		print '\tRetrieving results for chromosome number %s' % cnt
	for thrd in xrange(nThrds):
		qInJobs.put('STOP')
	#print errors
	if lErrsGlbl:
		print '\tThe following sequence(s) for gene(s) and position(s) were emptied:'
		print '\t%s' % '\n\t'.join(lErrsGlbl)
		print '\tAlignments with gaps in these positions were retrieved.'
	if lExcdGlbl:
		print '\tThe following interval(s) for gene(s) and position(s) were out of MAF:'
		print '\t%s' % '\n\t'.join(lExcdGlbl)
		print '\tAlignments for the interval(s) were not retrieved.'
	return 0


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest genome-wide in TargetMiner format.
########################################################
def procssUTRTrgtMnr(lDtObjcts, pthToMAF, extnsnMAF, extnsnUTRTrgtMnr, \
	pthFldrUTRFrmtdFlsTrgtMnr, sppRef, pthToUCSCKEGG, nThrds=5):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aMrdgIntrvls", "chr", and "strnd" for the genes 
	or lncrnas. pthToMAF is the path to the folder with the input MAF 
	files for each chromosome. extnsnMAF is the extension for the MAF 
	files in the pthToMAF path. pthFldrUTRFrmtdFlsTrgtMnr is the path to 
	the folder to output the gene/lncrna sequences in TargetMiner format 
	for each chromosome. sppRef is the UCSC name of the species (i.e. 
	hg19). nThrds is the number of parallel processes to run TrgtMnr. 
	extnsnUTRTrgtMnr is the extension of the output UTR files. 
	pthToUCSCKEGG is the path to the files relating all UCSC codes to
	other codes.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", and "seq" in the same order as the sorted names of 
	for the genes/lncrnas of interest in lDtObjcts. It also writes the
	TargetMiner formated sequence files in pthFldrUTRFrmtdFlsTrgtMnr. 
	Only those genes/lncrnas with information are going to be included 
	in the output. A list of errors in printed (if any).
	NOTE: A high number of threads requires a higher memory usage.
	"""
	#
	def wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFlsTrgtMnr, \
		extnsnUTRTrgtMnr, sppRef):
		""" Input: dGendMrgdNmSeqs is dictionary of {gene:{spp:seq}} 
		with the sequences merged, pthFldrUTRFrmtdFlsTrgtMnr is the 
		output folder to save the alignments, extnsnUTRTrgtMnr is
		the extension of the alignment file, sppRef is the reference 
		species. 
		This method will write the alignment sequences into the output
		folder in TargetMiner format.
		"""
		sGnNms = set(dGendMrgdNmSeqs.keys())#set of gene names
		while sGnNms:
			name = sGnNms.pop()
			outFlTrgtMnrFl = os.path.join(pthFldrUTRFrmtdFlsTrgtMnr, \
			''.join([name, extnsnUTRTrgtMnr]))
			dSppAlgndSeqs = dGendMrgdNmSeqs.pop(name)
			refSeq = lower(dSppAlgndSeqs.pop(sppRef)).replace('-', '')
			lFrtmdUTR = ['>%s %s\n%s\n>' % (name, name, refSeq)]
			#write file
			outFlTrgtMnrFl = open(outFlTrgtMnrFl, 'w')
			outFlTrgtMnrFl.write(''.join(lFrtmdUTR))
			outFlTrgtMnrFl.close()
		return 0
	#--------------------------
	#Start variables to test
	sUcscSpp = set([ l.split()[0] for l in open(pthToUCSCKEGG).read(). \
	splitlines() if l.strip() ])
	#Retrieve dictionary of sorted information by chromosome
	dChrlStrtEndGnNm, dGnNmStrnd = mkdChrlStrtEndGnNmStrnd(lDtObjcts, \
	pthFldrUTRFrmtdFlsTrgtMnr, extnsnUTRTrgtMnr)
	#Clear the info for chromosomes not present in pthToMAF
	sChrs = set(dChrlStrtEndGnNm.keys())#set of total chromosomes
	sChrsinMAFpth = set([ f.split(extnsnMAF)[0][3:] for f in \
	os.listdir(pthToMAF) if f.find(extnsnMAF) > -1 ])
	sChrAbsnt = sChrs.difference(sChrsinMAFpth)#real
	while sChrAbsnt:
		eChr = sChrAbsnt.pop()
		inf = dChrlStrtEndGnNm.pop(eChr)
		del inf
	#Retrieve alignments
	sChrsinMAFpth = sChrsinMAFpth.intersection(sChrs)
	lenChrs = len(sChrsinMAFpth)
	lErrsGlbl, lExcdGlbl = [], []
	qInJobs = Queue()
	qOutRslts = Queue()
	cnt = 0
	while sChrsinMAFpth:
		cnt += 1
		eChr = sChrsinMAFpth.pop()
		print ''.join(['\t...retrieving alignments for TargetMiner from', \
		' UTRs in chromosome %s, number %s out of %s' % (eChr, cnt, lenChrs)])
		fMAF = os.path.join(pthToMAF, 'chr%s.maf.gz' % eChr)
		lStrtEndGnNm = dChrlStrtEndGnNm.pop(eChr)
		qInJobs.put((fMAF, lStrtEndGnNm))
	for thrd in xrange(nThrds):
		Process(target=mthrdRtrvMAFs, args=(qInJobs,
		 qOutRslts,
		 sUcscSpp,
		 sppRef,
		 dGnNmStrnd,
		 pthFldrUTRFrmtdFlsTrgtMnr,
		 extnsnUTRTrgtMnr,
		 wrtAlgnmnt)).start()
	cnt = 0
	for eChr in xrange(lenChrs):
		cnt += 1
		lErrs, lExcd = qOutRslts.get()
		lErrsGlbl.extend(lErrs)
		lExcdGlbl.extend(lExcd)
		print '\tRetrieving results for chromosome number %s' % cnt
	for thrd in xrange(nThrds):
		qInJobs.put('STOP')
	#print errors
	if lErrsGlbl:
		print '\tThe following sequence(s) for gene(s) and position(s) were emptied:'
		print '\t%s' % '\n\t'.join(lErrsGlbl)
		print '\tAlignments with gaps in these positions were retrieved.'
	if lExcdGlbl:
		print '\tThe following interval(s) for gene(s) and position(s) were out of MAF:'
		print '\t%s' % '\n\t'.join(lExcdGlbl)
		print '\tAlignments for the interval(s) were not retrieved.'
	return 0


########################################################
#~ Make a dictionary of spp names as keys and sequences as values.
########################################################
def rtrvSubSeq(dNmSeqs, cStrtEnd, sppRef, trgtStrtEnd):
	""" Input: dNmSeqs is a dictionary of stored {names:sequences}, 
	cStrtEnd is a tuple of the genome start and end of the current 
	sequence of reference, sppRef is the species of reference, 
	trgtStrtEnd is the genome start and end of the sequence of interest. 
	Output: dSppAlgndSeqs is a dictionary of spp names as keys and 
	sequences as values. 
	Excluding only those that have only gaps. 
	NOTE: if reference only has gaps, then only the reference with 
	gaps will be returned"""
	cStrt, cEnd = cStrtEnd
	trgtStrt, trgtEnd = trgtStrtEnd
	strt = trgtStrt - cStrt
	end = trgtEnd - cStrt
	seq = dNmSeqs[sppRef]
	if set(seq) == {'-'}:#reference has only gaps
		seq = '-' * (end - strt)
		dSppAlgndSeqs = {sppRef: seq}
		return dSppAlgndSeqs
	else:
		srtdNONGapPos = [pos for pos in xrange(len(seq)) if seq[pos] \
		!= '-']
		#this might produce a bug but is added in order to consider the 
		#last nucleotide in the human alignment 
		srtdNONGapPos.append(pos + 1)
		strt = srtdNONGapPos[strt]
		end = srtdNONGapPos[end]
		dSppAlgndSeqs = {}
		for spp in dNmSeqs.keys():
			seq = dNmSeqs[spp][strt:end]
			if set(seq) != {'-'}:#not only gaps
				dSppAlgndSeqs[spp] = seq
		return dSppAlgndSeqs


########################################################
#~ Retrieve sequence alignments for genes/lncrnas or regions of 
# interest for one chromosome.
########################################################
def runCore(fMAF, lStrtEndGnNm, sUcscSpp, sppRef, dGnNmStrnd, \
	pthFldrUTRFrmtdFls, extnsnUTR, wrtAlgnmnt, TforU=False, \
	rtrndGendMrgdNmSeqs=False, dUcscTOMultizSpp=False):
	"""Input: fMAF is the full path to the MAF file. lStrtEndGnNm is a 
	list of start, end, and gene names. sUcscSpp is a set of species to 
	be included in the alignment (must include all in the MAF file). 
	sppRef is a reference species. dGnNmStrnd is a dictionary of 
	{gene name:strand}. pthFldrUTRFrmtdFls is the outfolder to save the 
	alignments, extnsnUTR is the extension of alignments output files. 
	wrtAlgnmnt is the method to write the alignment in a particular 
	format. Optionally, TforU return U instead of T, and if 
	rtrndGendMrgdNmSeqs a dictionary of gene names and sequences is 
	going to be returned instead of writing a file. Optionally, 
	dUcscTOMultizSpp is dictionary of multiz species and TargetScan 
	codes.
	Output: lErrs is a list for which part of the alignment for 
	the reference gene was found to be empty. lExcd is a list of 
	excluded intervals outside MAF files. Optionally, 
	rtrndGendMrgdNmSeqs a dictionary of gene names and sequences is 
	going to be returned instead of writing a file.
	"""
	#
	def overlap(start1, end1, start2, end2):
		return end1 >= start2 and end2 >= start1
	#open file
	if os.path.exists(fMAF):
		fMAF = gzip.open(fMAF, 'rb')
	else:
		fMAF = open(fMAF.split('.gz')[0], 'r')
	#Sort input list of positions
	lStrtEndGnNm.sort()
	#start variables
	blckToAppnd = []#make a list for the block holder
	strt, end, gnNm = lStrtEndGnNm.pop(0)
	trgtStrtEnd = (strt, end)
	cStrtEnd = (0, 0)
	dNmSeqs = dict([ (spp, '') for spp in sUcscSpp ])
	dGendIntrvldNmSeqs = {}#holder the exon alignments for each gene
	chckLstBlck = True
	#run core
	for elfMAF in fMAF:
		if not strt:#all sequences retrieved
			chckLstBlck = False
			break
		elif elfMAF.strip():
			if elfMAF[:7] == 'a score':
				dNmSeqs, cStrtEnd = apndSeqs(dNmSeqs, cStrtEnd, sppRef, \
				blckToAppnd, sUcscSpp)
				blckToAppnd = []
				cStrt, cEnd = cStrtEnd
				assert cStrt <= strt
				if overlap(strt, end, cStrt, cEnd):
					while end < cEnd:
						dSppAlgndSeqs = rtrvSubSeq(dNmSeqs, cStrtEnd, \
						sppRef, trgtStrtEnd)
						if dGendIntrvldNmSeqs.has_key(gnNm):
							dGendIntrvldNmSeqs[gnNm][strt, end] = \
							dSppAlgndSeqs
						else:
							dGendIntrvldNmSeqs[gnNm] = {(strt, end): \
							dSppAlgndSeqs}
						if lStrtEndGnNm:
							nxtStrt = lStrtEndGnNm[0][0]
							if nxtStrt > cEnd:#next start don't overlap 
								#with current interval
								dNmSeqs = dict([ (spp, '') for spp in \
								sUcscSpp ])
								cStrtEnd = (cEnd, cEnd)#reboot block
							strt, end, gnNm = lStrtEndGnNm.pop(0)
							trgtStrtEnd = (strt, end)
						else:#all sequences retrieved
							strt, end, gnNm = False, cEnd + 1, False
				else:#strt>cEnd
					dNmSeqs = dict([ (spp, '') for spp in sUcscSpp ])
					cStrtEnd = (cEnd, cEnd)#reboot block
			else:
				blckToAppnd.append(elfMAF)
	#include last block
	lExcd,sGnEcd = [],set()
	if chckLstBlck:
		dNmSeqs, cStrtEnd = apndSeqs(dNmSeqs, cStrtEnd, sppRef, \
		blckToAppnd, sUcscSpp)
		cStrt, cEnd = cStrtEnd
		while gnNm:
			if end >= cEnd:
				lExcd.append('\t'.join([gnNm, str(strt), str(end)]))
				sGnEcd.add(gnNm)
			else:
				dSppAlgndSeqs = rtrvSubSeq(dNmSeqs, cStrtEnd, sppRef, \
				trgtStrtEnd)
				if dGendIntrvldNmSeqs.has_key(gnNm):
					dGendIntrvldNmSeqs[gnNm][strt, end] = dSppAlgndSeqs
				else:
					dGendIntrvldNmSeqs[gnNm] = {(strt, end): \
					dSppAlgndSeqs}
			if lStrtEndGnNm:
				strt, end, gnNm = lStrtEndGnNm.pop(0)
				trgtStrtEnd = (strt, end)
			else:
				strt, end, gnNm = False, cEnd + 1, False
	#joined sequences
	dGendMrgdNmSeqs = {}
	sGns = set(dGendIntrvldNmSeqs.keys())
	lErrs = []#hold the list of genes and intervals that only have gaps
	while sGns:
		gn = sGns.pop()
		dIntrvldNmSeqs = dGendIntrvldNmSeqs.pop(gn)
		strnd = dGnNmStrnd[gn]
		dMrgdNmSeqs, lErrs = assmblGnExns(dIntrvldNmSeqs, sppRef, strnd, \
		gn, TforU)
		lErrs.extend(lErrs)
		dGendMrgdNmSeqs[gn] = dMrgdNmSeqs
	#exclude errors
	sGnEcd.update(set([x.split()[0] for x in lErrs]))
	for ecd in sGnEcd:
		if dGendMrgdNmSeqs.has_key(ecd):#in case part of the interval
			excd = dGendMrgdNmSeqs.pop(ecd)
	#return sequence
	if rtrndGendMrgdNmSeqs:
		return dGendMrgdNmSeqs, lErrs, lExcd
	else:
		if dUcscTOMultizSpp:
			wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFls, extnsnUTR, \
			sppRef, dUcscTOMultizSpp)
		else:
			wrtAlgnmnt(dGendMrgdNmSeqs, pthFldrUTRFrmtdFls, extnsnUTR, \
			sppRef)
		return lErrs, lExcd


########################################################
#~ Update a dictionary of spp names as keys and sequences as values.
########################################################
def updtSeqs(dNmSeqs, cStrtEnd, sppRef, trgtStrt):
	"""Input: dNmSeqs is a dictionary of stored {names:sequences}, 
	cStrtEnd tuple of the genome start and end of the current sequence 
	of reference, sppRef is the species of reference, trgtStrt is the 
	genome start and end of the sequence of interest.
	Output: dNmSeqs is the updated dictionary of stored 
	{names:sequences} starting on trgtStrt, cStrtEnd is the updated 
	genome start and end of the stored sequences.
	"""
	cStrt, end = cStrtEnd
	strt = trgtStrt - cStrt
	seq = dNmSeqs[sppRef]
	srtdNONGapPos = [pos for pos in xrange(len(seq)) if seq[pos] != '-']
   	#this might produce a bug but is added in order to consider the 
   	#last nucleotide in the human alignment 
	srtdNONGapPos.append(pos + 1)
	strt = srtdNONGapPos[strt]
	cStrtEnd = (strt, end)
	for spp in dNmSeqs.keys():
		seq = dNmSeqs[spp][strt:]
		dNmSeqs[spp] = seq
	return dNmSeqs, cStrtEnd

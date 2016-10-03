#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  runTest_alignments.py part of cellery (ceRNAs linking inference)
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
Test the alignment methods.
"""

__doc__ = "NOTE: multiprocessing must be patched as indicated by \
https://hg.python.org/cpython/rev/c82588ca3a79#l1.1"


########################################################
#~ Import libraries.
########################################################
from cellery import __object__
from cellery import alignments
from cellery import __classes__
from numpy import array,float32
from string import upper

import os
import unittest

#--------------------------------------------------------
#~ make dummy object
def mkLObjct_Isfrm():
	"""
	Make dummy object
	"""
	name = 'ENST00000536489'
	strnd = '-1'
	cmmNm = 'RPH3AL'
	pos = 0
	chr = '17'
	aIntrvls = array([(63435,63714),(65444,65593),(69413,69527),(96901, \
	97076),(169210,169340),(171062,171206),(177257,177370),(202502, \
	202633)])
	lenIntrvl = __object__.clcLen(aIntrvls)
	srtdMirnaNms = ['let-7a-5p','miR-181b-5p','miR-19b-2-5p']#sorted miRNAs
	aMirnaCnts = array([1,2,0],dtype=float32)#counts for TargetScan
	aMirnaRNAhCnts = array([1,2,1],dtype=float32)#counts for RNA hybrid
	aMirnaRNAhEngy = array([-21.3,-20.3,-20.0],dtype=float32) #energy for RNA hybrid
	aMirnaMirndCnts = array([3,2,0],dtype=float32)#counts for miRanda counts
	aMirnaMirndScrs = array([140.0,146.1,142.0],dtype=float32)#scores for miRanda 
	aMirnaMirndEngy = array([-21.3,-18.64,-15.64],dtype=float32)#energy for miRanda 
	aMirnaSVMicroCnts = array([1,2,3],dtype=float32)#counts for SVMicro
	aMirnaSVMicroScrs = array([240.0,136.1,122.5],dtype=float32)#scores for miRanda 
	aMirnaTrgtMnrCnts = array([1,0,3],dtype=float32)#counts for TargetMiner
	aMirnaPITACnts = array([1,1,3],dtype=float32)#counts for miRanda 
	aMirnaPITAScrs = array([40.0,246.1,42.0],dtype=float32)#scores for miRanda 
	aMirnaMirMapScrs = array([32.0,44.1,21.0],dtype=float32)#scores for mirMap
	aMirnaMirMapCnts = array([2,2,3],dtype=float32)#counts for mirMap
	aMirnaWalk = array([0,1,0],dtype=bool)#mirWalk confirmation
	#make gen object
	gene = __classes__.gene(name,srtdMirnaNms)
	gene.strnd = strnd
	gene.cmmNm = cmmNm
	gene.pos = pos
	gene.chr = chr
	gene.aIntrvls = aIntrvls
	gene.len = lenIntrvl
	gene.aMirnaCnts = aMirnaCnts
	gene.aMirnaRNAhCnts = aMirnaRNAhCnts
	gene.aMirnaRNAhEngy = aMirnaRNAhEngy
	gene.aMirnaMirndCnts = aMirnaMirndCnts
	gene.aMirnaMirndScrs = aMirnaMirndScrs
	gene.aMirnaMirndEngy = aMirnaMirndEngy
	gene.aMirnaSVMicroCnts = aMirnaSVMicroCnts
	gene.aMirnaSVMicroScrs = aMirnaSVMicroScrs
	gene.aMirnaTrgtMnrCnts = aMirnaTrgtMnrCnts
	gene.aMirnaPITACnts = aMirnaPITACnts
	gene.aMirnaPITAScrs = aMirnaPITAScrs
	gene.aMirnaMirMapScrs = aMirnaMirMapScrs
	gene.aMirnaMirMapCnts = aMirnaMirMapCnts
	gene.aMirnaWalk = aMirnaWalk
	return [gene]


########################################################
#~ Retrieve sequence alignments in fasta/mirMap format.
########################################################
def runTest_fastaMirMap():
	#--------------------------------------------------------
	#~ Define variables
	dbFldr = os.path.join('tests2','db')# db folder
	pthToMAF = os.path.join(dbFldr,'maf')# multiz100-way MAF files
	extnsnMAF = '.maf' # extension for the MAF files
	extnsnUTRMirMap = '.fas'# extension for the input alignments for mirmap
	sppRef = 'hg19' #spp code for UCSC alignments (i.e. a name==hg19)
	pthToUCSCKEGG = os.path.join(dbFldr,'dbUCSCtKEGG.tsv')
	lObjcts_Isfrm = mkLObjct_Isfrm()
	#--------------------------------------------------------
	#~ Define output folders
	alignmntsFldr = os.path.join(dbFldr,'alignments')
	if not os.path.exists(alignmntsFldr):
		os.mkdir(alignmntsFldr)
	#--------------------------------------------------------
	#~ Return alignments for protein coding isoforms
	alignments.proccsMAFnUTRsMirMap(lObjcts_Isfrm,pthToMAF,extnsnMAF, \
	extnsnUTRMirMap,alignmntsFldr,sppRef,pthToUCSCKEGG)
	#--------------------------------------------------------
	#~ Test result
	inFst = os.path.join(alignmntsFldr,'%s%s'%(lObjcts_Isfrm[0].name, \
	extnsnUTRMirMap))
	seqTest = dict([(seq.split()) for seq in open(inFst).read(). \
	split('>') if seq.strip()])[sppRef]
	seqTest = seqTest.replace('-','')
	return seqTest
	

########################################################
#~ Retrieve sequence alignments in fasta/RNAhybrid/PITA/miRanda format 
# (fasta only reference).
########################################################
def runTest_fastaRNAhybridPITAmiRanda():
	#--------------------------------------------------------
	#~ Define variables
	dbFldr = os.path.join('tests2','db')# db folder
	pthToMAF = os.path.join(dbFldr,'maf')# multiz100-way MAF files
	extnsnMAF = '.maf' # extension for the MAF files
	extnsnUTRMirMap = '.fasta'# extension for the input alignments for mirmap
	sppRef = 'hg19' #spp code for UCSC alignments (i.e. a name==hg19)
	pthToUCSCKEGG = os.path.join(dbFldr,'dbUCSCtKEGG.tsv')
	lObjcts_Isfrm = mkLObjct_Isfrm()
	#--------------------------------------------------------
	#~ Define output folders
	alignmntsFldr = os.path.join(dbFldr,'alignments')
	if not os.path.exists(alignmntsFldr):
		os.mkdir(alignmntsFldr)
	#--------------------------------------------------------
	#~ Return alignments for protein coding isoforms
	alignments.proccsMAFnUTRsOnlySppRef(lObjcts_Isfrm,pthToMAF, \
	extnsnMAF,extnsnUTRMirMap,alignmntsFldr,sppRef,pthToUCSCKEGG)
	#--------------------------------------------------------
	#~ Test result
	inFst = os.path.join(alignmntsFldr,'%s%s'%(lObjcts_Isfrm[0].name, \
	extnsnUTRMirMap))
	seqTest = dict([(seq.split()) for seq in open(inFst).read(). \
	split('>') if seq.strip()])['ENST00000536489']
	seqTest = seqTest.replace('-','')
	return seqTest


########################################################
#~ Retrieve sequence alignments in TargetMiner format.
########################################################
def runTest_targetMiner():
	#--------------------------------------------------------
	#~ Define variables
	dbFldr = os.path.join('tests2','db')# db folder
	pthToMAF = os.path.join(dbFldr,'maf')# multiz100-way MAF files
	extnsnMAF = '.maf' # extension for the MAF files
	extnsnUTRTrgtMnr = '.tmr'# extension of UTR file for TargetMiner
	sppRef = 'hg19' #spp code for UCSC alignments (i.e. a name==hg19)
	pthToUCSCKEGG = os.path.join(dbFldr,'dbUCSCtKEGG.tsv')
	lObjcts_Isfrm = mkLObjct_Isfrm()
	#--------------------------------------------------------
	#~ Define output folders
	alignmntsFldr = os.path.join(dbFldr,'alignments')
	if not os.path.exists(alignmntsFldr):
		os.mkdir(alignmntsFldr)
	#--------------------------------------------------------
	#~ Return alignments for protein coding isoforms
	alignments.procssUTRTrgtMnr(lObjcts_Isfrm,pthToMAF,extnsnMAF, \
	extnsnUTRTrgtMnr,alignmntsFldr,sppRef,pthToUCSCKEGG)
	#--------------------------------------------------------
	#~ Test result
	inFst = os.path.join(alignmntsFldr,'%s%s'%(lObjcts_Isfrm[0].name, \
	extnsnUTRTrgtMnr))
	seqTest = dict([(seq.split()[0],seq.splitlines()[1]) for seq in \
	open(inFst).read().split('>') if seq.strip()])['ENST00000536489']
	seqTest = upper(seqTest.replace('-',''))
	return seqTest


########################################################
#~ Retrieve sequence alignments in TargetScan format.
########################################################
def runTest_targetScan():
	#--------------------------------------------------------
	#~ Define variables
	dbFldr = os.path.join('tests2','db')# db folder
	pthToMAF = os.path.join(dbFldr,'maf')# multiz100-way MAF files
	extnsnMAF = '.maf' # extension for the MAF files
	extnsnUTRTrgtScn = '.utr'# extension of UTR file for TargetScan
	sppRef = 'hg19' #spp code for UCSC alignments (i.e. a name==hg19)
	pthToUCSCKEGG = os.path.join(dbFldr,'dbUCSCtKEGG.tsv')
	lObjcts_Isfrm = mkLObjct_Isfrm()
	#--------------------------------------------------------
	#~ Define output folders
	alignmntsFldr = os.path.join(dbFldr,'alignments')
	if not os.path.exists(alignmntsFldr):
		os.mkdir(alignmntsFldr)
	#--------------------------------------------------------
	#~ Return alignments for protein coding isoforms
	alignments.proccsMAFnUTRsTrgtScn(lObjcts_Isfrm,pthToMAF,extnsnMAF, \
	extnsnUTRTrgtScn,alignmntsFldr,sppRef,pthToUCSCKEGG)
	#--------------------------------------------------------
	#~ Test result
	inFst = os.path.join(alignmntsFldr,'%s%s'%(lObjcts_Isfrm[0].name, \
	extnsnUTRTrgtScn))
	seqTest = dict([(seq.split()[1],seq.split()[2]) for seq in \
	open(inFst).read().split('>') if seq.strip()])['9606']
	seqTest = seqTest.replace('-','')
	seqTest = seqTest.replace('U','T')
	return seqTest


########################################################
#~ Retrieve sequence alignments in SVMicro format.
########################################################
def runTest_SVMicro():
	#--------------------------------------------------------
	#~ Define variables
	dbFldr = os.path.join('tests2','db')# db folder
	pthToPhstCns = os.path.join(dbFldr,'phastCons46way.d')
	extnsnPhstCns = '.phastCons46way.wigFix.gz'
	pthToMAF = os.path.join(dbFldr,'maf')# multiz100-way MAF files
	extnsnMAF = '.maf' # extension for the MAF files
	extnsnSVMicro = '.svm'# extension of UTR file for SVMicro
	sppRef = 'hg19' #spp code for UCSC alignments (i.e. a name==hg19)
	pthToUCSCKEGG = os.path.join(dbFldr,'dbUCSCtKEGG.tsv')
	lObjcts_Isfrm = mkLObjct_Isfrm()
	#--------------------------------------------------------
	#~ Define output folders
	alignmntsFldr = os.path.join(dbFldr,'alignments')
	if not os.path.exists(alignmntsFldr):
		os.mkdir(alignmntsFldr)
	#--------------------------------------------------------
	#~ Return alignments for protein coding isoforms
	alignments.procssUTRSVMicro(lObjcts_Isfrm,pthToPhstCns, \
	extnsnPhstCns,alignmntsFldr,extnsnSVMicro,pthToMAF,extnsnMAF, \
	sppRef,pthToUCSCKEGG)
	#--------------------------------------------------------
	#~ Test result
	inFst = os.path.join(alignmntsFldr,'%s%s'%(lObjcts_Isfrm[0].name, \
	extnsnSVMicro))
	seqTest = dict([(seq.splitlines()[1].split()[1], \
	seq.splitlines()[1].split()[3]) for seq in open(inFst).read(). \
	split('>') if seq.strip()])['ENST00000536489']
	return seqTest


#--------------------------------------------------------
#~ run test
class TestAlignmentMethods(unittest.TestCase):
	global trueSeq
	trueSeq = ''.join(['ATTGACTTAAGTCCCAGTGATTCAGCTCCTCATCTGGAA', \
	'CACCTCGGGTCACCCCCGACAACGGTGGTGGGAGGGAGAGCGGCCTCCTCCTCCCTGGTGGGG', \
	'CCTGTCTGGGTGAAGCCCCTCTGTTCCCGATGTGACTCCCCACCCCCAGCCGGGTGCTCCGAG', \
	'CCATGGCCGACACCATCTTCGGCAGCGGGAATGATCAGTGGGTTTGCCCCAATGACCGGCAGC', \
	'TTGCCCTTCGAGCCAAGCTGCAGACGGGCTGGTCCGTGCACACCTACCAGACGGAGAAGCAGA', \
	'GGAGGAAGCAGCACCTCAGCCCGGCGGAGGTGGAGGCCATCCTGCAGGTCATCCAGAGGGCAG', \
	'AGCGGCTCGACGTCCTGGAGCAGCAGAGAATCGGGCGGCTGGTGGAGCGGCTGGAGACCATGA', \
	'GGCGGAATGTGATGGGGAACGGCCTGTCCCAGTGTCTGCTCTGCGGGGAGGTGCTGGGCTTCC', \
	'TGGGCAGCTCGTCGGTGTTCTGCAAAGACTGCAGGAAGGTCTGGAAGAGGTCGGGGGCCTGGT', \
	'TCTACAAAGGGCTCCCCAAGTATATCTTGCCCCTGAAGACCCCTGGCCGAGCTGATGACCCCC', \
	'ACTTCCGACCTTTGCCCACGGAACCGGCAGAGCGAGAGCCCAGAAGCTCTGAGACCAGCCGCA', \
	'TCTACACGTGGGCCCGAGGAAGAGTGGTTTCCAGTGACAGTGACAGTGACTCGGATCTTAGCT', \
	'CCTCCAGCCTAGAGGACAGACTCCCATCCACTGGGGTCAGGGACCGGAAAGGCGACAAACCCT', \
	'GGAAGGAGTCAGGTGGCAGCGTGGAGGCCCCCAGGATGGGGTTCACCCACCCGCCGGGCCACC', \
	'TCTCTGGGTGCCAGAGCAGCCTGGCCAGTGGTGAGACGGGGACAGGCTCTGCTGACCCGCCAG', \
	'GGGGACCCCGCCCCGGGCTGACCCGAAGGGCCCCGGTAAAAGACACACCTGGACGAGCCCCCG', \
	'CTGCTGACGCAGCTCCAGCAGGCCCCTCCAGCTGCCTGGGCTGAGGTGTCTGGTGCCTGGAAC', \
	'AGACTTCCCTGTGGAGGATTCCTGCCAGACCCTGCCCGGCTCCTCCCTGACCGGTCCTTGTGC', \
	'CCTCACCAGACACCCTGTTGGCCATGACTCAACAAACCAGTGTTGGGAGCCGTCTGCCTCCCC', \
	'AGCTCAGTGCCTTTCTGCACCCCTTCTCTCCTGGGGAGCTGTCTGCATCCGCCACCCCCTCC'])
	#--------------------------------------------------------
	def test_fastaMirMap(self):
		self.assertEqual(runTest_fastaMirMap(),trueSeq)
	def test_fastaRNAhybridPITAmiRanda(self):
		self.assertEqual(runTest_fastaRNAhybridPITAmiRanda(),trueSeq)
	def test_targetMiner(self):
		self.assertEqual(runTest_targetMiner(),trueSeq)
	def test_targetScan(self):
		self.assertEqual(runTest_targetScan(),trueSeq)
	def test_SVMicro(self):
		self.assertEqual(runTest_SVMicro(),trueSeq)
	
		
if __name__ == '__main__':
    unittest.main()


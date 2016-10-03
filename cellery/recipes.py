#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  recipes.py part of cellery (ceRNAs linking inference)
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
Common recipes for analysis
"""

########################################################
#~ Import libraries.
########################################################
from cellery import __object__
from cellery import coexpression
from cellery import ecp
from cellery import exceptions
from cellery import matrices
from cellery import networks
from cellery import randomz
from cellery import sampling
from cellery import statistics
from cellery import wrapper
from numpy import array,float32,ma,newaxis,ones

#----------------------------
#Recipes to format input
#----------------------------

########################################################
#~ Construct objects for sequences of protein coding genes/isoforms, 
# 3'UTRs genes/isoforms, and lncRNAs genes/isoforms from sql databases.
########################################################
def mklObjcts3UTRLncrns(fOutPCSql,fOutUtr3cdgSql,fOutLncrnaSql, \
	srtdMirnaFms):
	"""
	Input: fOutPCSql is a sql database with information for protein 
	coding genes/isoforms objects about "pos", "name", "aIntrvls", 
	"chr", "cmmNm", "strnd". fOutUtr3cdgSql is a sql database with 
	information for 3'UTRs of protein coding genes/isoforms objects 
	about "pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd". 
	fOutLncrnaSql is a sql database with information for lncRNA genes/
	isoforms objects about "pos", "name", "aIntrvls", "chr", "cmmNm", 
	and "strnd". srtdMirnaFms is the sorted miRNAs of interest.
	Output: lObjctsPC is a list of objects for protein coding  genes/
	isoforms with the pointers "pos", "name", "aIntrvls", "chr", 
	"cmmNm", and "strnd". lObjctsUtr3cdg  is a list of objects for 
	3'UTRS of protein coding genes/isoforms with the pointers "pos", 
	"name", "aIntrvls", "chr", "cmmNm", and "strnd". lObjctsLncrna is a 
	list of objects for lncRNA coding genes/isoforms with the pointers 
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd".
	"""
	#--------------------------------------------------------
	#~ Build objects for protein coding genes/isoforms
	lObjctsPC = __object__.mkArrySqlObjcts(fOutPCSql,srtdMirnaFms)
	#~ Build objects for all 3'UTRs from protein coding genes/isoforms
	lObjctsUtr3cdg = __object__.mkArrySqlObjcts(fOutUtr3cdgSql, \
	srtdMirnaFms)
	#~ Build objects for a subset of lncRNA genes/isoforms -All 
	# annotated by GENCODE-
	lObjctsLncrna = __object__.mkArrySqlObjcts(fOutLncrnaSql, \
	srtdMirnaFms)
	return lObjctsPC,lObjctsUtr3cdg,lObjctsLncrna


#----------------------------
#Recipes to run and build object with MREs from all prediction programs
#----------------------------
########################################################
#~ Run all MRE prediction algoirthms for all the obecjts in a list and
# update the predicted values
########################################################
"""
extnsnMAF is the extension for the MAF files in the pthToMAF path. 
extnsnMirMap is the extension of mirMap file results. extnsnMiranda is 
the extension of miRanda file results. extnsnMirnTrgtMnr is the 
extension of the mirna TargetMiner-formated output file. 
extnsnMirnaTrgtScn is the extension of the mirna TargetScan-formated 
(fasta) output file. extnsnPITA is the extension of PITA file results.
extnsnPhstCns is the extension for the phastCons score files in the 
pthToPhstCns path. extnsnPrwsTrgtMnr is the extension of the TargetMiner-
formated pair files. extnsnRNAhybrd is the extension of RNAhybrid file 
results. extnsnSVMicro is the extension of the output UTR files (in 
SVMicro format). extnsnTrgtMnr is the extension to the output files from 
TargetMiner. extnsnTrgtScn is the extension of TargetScan file results.
extnsnMAF is the extension for the MAF files in the pthToMAF path. 
extnsnUTRTrgtMnr is the extension of the output UTR files.
extnsnUTRTrgtScn is the extension of the output UTR files. fOutDtSql is 
a sqlite3 database with the attributes of interest (for objects).
fstaDtFl is the name of the input fasta file (for objects).
fstasMirnaBsFl is an array of miRNA fasta-formated file names.mirbaseSpp 
is the species code in miRBase base (i.e. ==KEGG format: {hsa,mmu}).
mirnFrmtdFlTrgtMnr is the output file to write the TargetMiner-formatted 
output file. nPrlProcss is the number of parallel processes to run.
nThrds is the number of cores to run. phastModFl is the phastCons model 
file (can be obtained from the UCSC browser multiZ version of interest).
phastTree is the phastCons tree file (can be obtained from the UCSC 
browser multiZ version of interest). pthFldrFstaDtFl is the path to the 
folder with fasta sequences for objects. pthFldrFstaMirnsFls is a path 
to the folder with the sequence files of the objects. 
pthFldrMirnsFrmtdFlsTrgtMnr is a path to the folder with the files with 
prefix mirnaFrmtdFlTrgtMnr. pthFldrMirnsFrmtdFlsTrgtScn is a path to the 
folder with the fasta files from miRNAs. pthFldrPrwsFrmtdFlsTrgtMnr is a 
path to the folder with the files in TargetMiner-formated pair files.
pthFldrUTRFrmtdFlsFsta is the path to the folder tooutput thegene/lncrna 
sequences in miRMap/fasta format for each chromosome. 
pthFldrUTRFrmtdFlsFstaOnlySppRef is the path to the folder to output the 
gene/lncrna sequences in OnlySppRef/fasta format for each chromosome.
pthFldrUTRFrmtdFlsSVMicro is the path to the folder to output the 
gene/lncrna sequences in SVMicro format. pthFldrUTRFrmtdFlsTrgtMnr isthe 
path to the folder to output the gene/lncrna sequences in TargetMiner 
format for each chromosome. pthFldrUTRFrmtdFlsTrgtScn is the path to the 
folder to output the gene/lncrna sequences in  fomra for each chromosome.
pthMirMapLib is the full path to the mirMap python library.pthMiranda is 
the path to the miranda executable. pthOutFldrMirMap is a path to the 
folder to save the MirMap results. pthOutFldrMiranda is a path to the 
folder to save the miRanda results. pthOutFldrPITA isa path tothe folder 
to save the PITA results. pthOutFldrRNAhybrd is a path to the folder to 
save the RNAhybrid results. pthOutFldrSVMicro is a path to the folder 
with the SVMicro results. pthOutFldrTrgtMnr is a path to the folder with 
the TargetMiner results. pthOutFldrTrgtScn is a path to the folder with 
the TargetScan results. pthPITA is the path the PITA executable.pthPhast 
is the full path to the phast executable. pthPython is the full path to 
python. pthPythonLib is the full path to the python library. pthRNAhybrd 
is the path to the RNAhybrid executable. pthSVMicro is the path to the 
SVMicro executable. pthTmpFldr is a folder to save temporal  files. 
pthToMAF is the path to the folder with the input MAF files for each 
chromosome. pthToPhstCns is the path to phastCons files.pthToUCSCKEGG is 
the path to the files relating all UCSC codes to other codes. pthTrgtMnr 
is the to the TargetMiner executable. pthTrgtScn is the path to the 
TargetScan executable. sppMltzCd is the species code of interest for 
MultiZ alignment (is a number).sppRef is a reference species.sqlFlMirMap 
sqlFl is a sqlite3 database with fields having (a subset) of names inthe 
same order of aMirnaNms for mirMap results.sqlFlMirnd sqlFl is a sqlite3 
database with fields having (a subset) of names in the same order of 
aMirnaNms for miRanda results.sqlFlPITA sqlFl is a sqlite3 database with 
fields having (a subset) of names in the same order of aMirnaNms forPITA 
results.sqlFlRNAhybrid sqlFl is a sqlite3 database with fields having (a 
subset) of names in the same order of aMirnaNms for RNAhybrid results.
sqlFlSVMicro sqlFl is a sqlite3 database with fields having (asubset) of 
names in the same order of aMirnaNms for SVMicro results. sqlFlTrgtMnr 
sqlFl is a sqlite3 database with fields having (a subset) of names inthe 
same order of aMirnaNms for TargetMiner results. sqlFlTrgtScn sqlFl is a 
sqlite3 database with fields having (a subset) of names inthe same order 
of aMirnaNms for TargetScaner results.
NOTE: Variables are set to human by default
"""
def mklObjctsRnAllPrdctns(pthFldrFstaDtFl, \#folder to the fasta file for the object of interest
	fstaDtFl, \#name of the fasta file for the object of interest
	rsltsFldr, \#results folder
	pthFldrFstaMirnsFls = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/GTeX_Mele.d/data.d/mirnas_data.d/hg_19.d', \#folders for miRNAs files for TrgtScn in fasta
	fstasMirnaBsFl = 'mature_hg19.fa', \#to the fasta file with the miRNA data in fasta format	
	pthTmpFldr = 'tmp.d', \#folder to store temporal files
	extnsnMAF = '.maf', \# extension for the MAF files
	extnsnMirMap = '.mmp', \# extension of output result files for mirMap
	extnsnMiranda = '.mrd', \# extension for the input sequences for miRanda
	extnsnMirnTrgtMnr = '.mmt', \# extension of miRNA file for TargetMiner
	extnsnMirnaTrgtScn = '.its', \# extension for miRNA targetscan files
	extnsnPITA = '.pta', \# extension of output result files for PITA
	extnsnPhstCns = '.phastCons46way.wigFix.gz', \#multiZ 46-way
	extnsnPrwsTrgtMnr = '.tpr', \# extension of pairwise comparisons file for TargetMiner
	extnsnRNAhybrd = '.rhy', \# extension for the input sequences for RNAhybrid
	extnsnSVMicro = '.svm', \# extension of output result files for SVMicro
	extnsnTrgtMnr = '.tgm', \# extension of output result files for TargetMiner
	extnsnTrgtScn = '.tsc', \# extension for targetscan files
	extnsnUTRFsta = '.fas', \# extension for the input alignments for fasta
	extnsnUTRTrgtMnr = '.tmr', \# extension of UTR file for TargetMiner
	extnsnUTRTrgtScn = '.utr', \# extension of UTR file for TargetScan
	mirbaseSpp = 'hsa', \#spp code for miRBase base (i.e. ==KEGG format: {hsa,mmu}). 
	nPrlProcss = 500, \#number of parallel processes to run in queue	
	nThrds = 5, \# number of threads for multiprocessing
	phastModFl = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/data/complex1_4.d/hg19.d/maf/hg19.100way.phastCons.mod', \
	phastTree = '((((((((((((((((((hg19,panTro4),gorGor3),ponAbe2),nomLeu3),(((rheMac3,macFas5),papHam1),chlSab1)),(calJac3,saiBol1)),otoGar3),tupChi1),(((speTri2,(jacJac1,((micOch1,(criGri1,mesAur1)),(mm10,rn5)))),(hetGla2,(cavPor3,(chiLan1,octDeg1)))),(oryCun2,ochPri3))),(((susScr3,((vicPac2,camFer1),((turTru2,orcOrc1),(panHod1,(bosTau7,(oviAri3,capHir1)))))),(((equCab2,cerSim1),(felCat5,(canFam3,(musFur1,(ailMel1,(odoRosDiv1,lepWed1)))))),((pteAle1,pteVam1),((myoDav1,myoLuc2),eptFus1)))),(eriEur2,(sorAra2,conCri1)))),(((((loxAfr3,eleEdw1),triMan1),(chrAsi1,echTel2)),oryAfe1),dasNov3)),(monDom5,(sarHar1,macEug2))),ornAna1),(((((((falChe1,falPer1),(((ficAlb2,((zonAlb1,geoFor1),taeGut2)),pseHum1),(melUnd1,(amaVit1,araMac1)))),colLiv1),(anaPla1,(galGal4,melGal1))),allMis1),((cheMyd1,chrPic1),(pelSin1,apaSpi1))),anoCar2)),xenTro7),latCha1),(((((((tetNig2,(fr3,takFla1)),(oreNil2,(neoBri1,(hapBur1,(mayZeb1,punNye1))))),(oryLat2,xipMac1)),gasAcu1),gadMor1),(danRer7,astMex1)),lepOcu1)),petMar2)', \
	pthMirMapLib = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/code.d/miRmap-1.1/libs/lib-archlinux-x86_64', \
	pthMiranda = '/exports/igmm/datastore/ponting-lab/software.d/bin', \
	pthPITA = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/code.d/PITA', \
	pthPhast = '/exports/igmm/datastore/ponting-lab/software.d/bin', \
	pthPython = '/exports/igmm/software/pkg/el7/apps/python/2.7.10/bin', \
	pthPythonLib='/exports/igmm/datastore/ponting-lab/software.d', \
	pthRNAhybrd = '/exports/igmm/datastore/ponting-lab/software.d/bin', \
	pthSVMicro = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/code.d/SVMicrO', \
	pthToMAF = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/data/complex1_4.d/hg19.d/maf', \
	pthToPhstCns = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/data/db.d/phastCons46way.d', \
	pthToUCSCKEGG = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/data/dbUCSCtKEGG.tsv', \
	pthTrgtMnr = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/code.d/TargetMiner-executable/TargetMiner', \
	pthTrgtScn = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/code.d/targetscan_70', \
	sppMltzCd = '9606', \#spp code for multiz alignments (a number)
	sppRef = 'hg19')#spp code for UCSC alignments (i.e. a name==hg19)
	#----------------------------
	# Define output files, databases and folders
	mirnFrmtdFlTrgtMnr = '%s_TrgtMnr'%os.path.splitext(fstasMirnaBsFl)
	alignmntsFldr = os.path.join(rsltsFldr,'algmnts.d')
	pthFldrUTRFrmtdFlsFsta = os.path.join(alignmntsFldr,'fasta.d')
	pthFldrUTRFrmtdFlsFstaOnlySppRef = os.path.join(alignmntsFldr, \
	'fastaOnlySppRef.d')
	pthFldrUTRFrmtdFlsSVMicro = os.path.join(alignmntsFldr,'SVMicro.d')
	pthFldrUTRFrmtdFlsTrgtMnr = os.path.join(alignmntsFldr,'trgtMnr.d')
	pthFldrUTRFrmtdFlsTrgtScn = os.path.join(alignmntsFldr,'trgtScn.d')
	pthFldrMirnsFrmtdFlsTrgtMnr =  os.path.join(rsltsFldr, \
	'trgtMnrMirns.d')
	pthFldrMirnsFrmtdFlsTrgtScn =  os.path.join(rsltsFldr, \
	'trgtScnMirns.d')
	pthFldrPrwsFrmtdFlsTrgtMnr =  os.path.join(rsltsFldr, \
	'trgtMnrPrws.d')
	pthOutFldrMirMap = os.path.join(rsltsFldr,'mirMap.d')
	pthOutFldrMiranda = os.path.join(rsltsFldr,'miRanda.d')
	pthOutFldrPITA = os.path.join(rsltsFldr,'PITA.d')
	pthOutFldrRNAhybrd = os.path.join(rsltsFldr,'RNAhybrd.d')
	pthOutFldrSVMicro = os.path.join(rsltsFldr,'SVMicro.d')
	pthOutFldrTrgtMnr = os.path.join(rsltsFldr,'trgtMnr.d')
	pthOutFldrTrgtScn = os.path.join(rsltsFldr,'trgtScn.d')
	fOutDtSql = os.path.join(rsltsFldr,'%s.db'%os.path. \
	splitext(fstaDtFl)[0])
	sqlFlMirMap = os.path.join(rsltsFldr,'mirMap.db')
	sqlFlMirnd = os.path.join(rsltsFldr,'miRanda.db')
	sqlFlPITA = os.path.join(rsltsFldr,'PITA.db')
	sqlFlRNAhybrid = os.path.join(rsltsFldr,'RNAhybrd.db')
	sqlFlSVMicro = os.path.join(rsltsFldr,'SVMicro.db')
	sqlFlTrgtMnr = os.path.join(rsltsFldr,'trgtMnr.db')
	sqlFlTrgtScn = os.path.join(rsltsFldr,'trgtScn.db')
	#----------------------------
	# Make folders
	for fldr in [alignmntsFldr,pthTmpFldr,rsltsFldr, \
	pthFldrUTRFrmtdFlsFsta,pthFldrUTRFrmtdFlsFstaOnlySppRef, \
	pthFldrUTRFrmtdFlsSVMicro,pthFldrUTRFrmtdFlsTrgtMnr, \
	pthFldrUTRFrmtdFlsTrgtScn,pthFldrMirnsFrmtdFlsTrgtMnr, \
	pthFldrMirnsFrmtdFlsTrgtScn,pthFldrPrwsFrmtdFlsTrgtMnr, \
	pthOutFldrMirMap,pthOutFldrMiranda,pthOutFldrPITA, \
	pthOutFldrRNAhybrd,pthOutFldrSVMicro,pthOutFldrTrgtMnr, \
	pthOutFldrTrgtScn]:
	if not os.system.exists(fldr):
		os.mkdir(fldr)	
	#----------------------------
	#Run TargetScan from fasta and MAF files / sql databases
	lObjcts_trgt = wrapper.runTrgtScn(extnsnMAF,extnsnMirnaTrgtScn, \
	extnsnTrgtScn,extnsnUTRTrgtScn,fstaDtFl,fstasMirnaBsFl,mirbaseSpp, \
	nPrlProcss,nThrds,pthFldrFstaDtFl,pthFldrFstaMirnsFls, \
	pthFldrMirnsFrmtdFlsTrgtScn,pthFldrUTRFrmtdFlsTrgtScn, \
	pthOutFldrTrgtScn,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	pthTrgtScn,sppRef,sppMltzCd,fOutDtSql,sqlFlTrgtScn)
	#----------------------------
	#Run miRanda from fasta and MAF files / sql databases and update 
	lDtObjcts_toAdd = wrapper.runMiRanda(extnsnMAF,extnsnMiranda, \
	extnsnUTRFsta,fstaDtFl,fstasMirnaBsFl,nPrlProcss,nThrds, \
	pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFstaOnlySppRef, \
	pthOutFldrMiranda,pthMiranda,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	sppRef,fOutDtSql,sqlFlMirnd)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaMirndCnts','aMirnaMirndScrs','aMirnaMirndEngy'])
	#----------------------------
	#Run mirMap from fasta and MAF files / sql databases  and update 
	lDtObjcts_toAdd = wrapper.runMirMap(extnsnMAF,extnsnMirMap, \
	extnsnUTRFsta,fstaDtFl,fstasMirnaBsFl,nPrlProcss,nThrds,phastModFl, \
	phastTree,pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFsta, \
	pthMirMapLib,pthOutFldrMirMap,pthPython,pthPythonLib,pthPhast, \
	pthTmpFldr,pthToMAF,pthToUCSCKEGG,sppRef,fOutDtSql,sqlFlMirMap)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaMirMapScrs','aMirnaMirMapCnts'])
	#----------------------------
	#Run PITA from fasta and MAF files / sql databases and update 
	lDtObjcts_toAdd = wrapper.runPITA(extnsnMAF,extnsnUTRFsta, \
	extnsnPITA,fstaDtFl,fstasMirnaBsFl,nPrlProcss,nThrds, \
	pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFstaOnlySppRef, \
	pthOutFldrPITA,pthPITA,pthTmpFldr,pthToMAF,pthToUCSCKEGG,sppRef, \
	fOutDtSql,sqlFlPITA)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaPITACnts','aMirnaPITAScrs'])
	#----------------------------
	#Run RNAhybrid from fasta and MAF files / sql databases and update 
	lDtObjcts_toAdd = wrapper.runRNAhybrid(extnsnMAF,extnsnRNAhybrd, \
	extnsnUTRFsta,fstaDtFl,fstasMirnaBsFl,nPrlProcss,nThrds, \
	pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsFstaOnlySppRef, \
	pthOutFldrRNAhybrd,pthRNAhybrd,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	sppRef,fOutDtSql,sqlFlRNAhybrid)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaRNAhCnts','aMirnaRNAhEngy'])
	#----------------------------
	#Run TargetMiner from fasta and MAF files / sql databases and update 
	lDtObjcts_toAdd = wrapper.runTargetMiner(extnsnMAF,extnsnMirnTrgtMnr, \
	extnsnPrwsTrgtMnr,extnsnTrgtMnr,extnsnUTRTrgtMnr,fstaDtFl, \
	fstasMirnaBsFl,mirnFrmtdFlTrgtMnr,nPrlProcss,nThrds,pthFldrFstaDtFl, \
	pthFldrFstaMirnsFls,pthFldrMirnsFrmtdFlsTrgtMnr, \
	pthFldrPrwsFrmtdFlsTrgtMnr,pthFldrUTRFrmtdFlsTrgtMnr, \
	pthOutFldrTrgtMnr,pthTrgtMnr,pthTmpFldr,pthToMAF,pthToUCSCKEGG, \
	sppRef,fOutDtSql,sqlFlTrgtMnr,lAttrbts)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaTrgtMnrCnts'])
	#----------------------------
	#Run SVMicro from fasta and MAF files / sql databases and update 
	lDtObjcts_toAdd = wrapper.runSVMicro(extnsnMAF,extnsnPhstCns, \
	extnsnSVMicro,fstaDtFl,fstasMirnaBsFl,nPrlProcss,nThrds, \
	pthFldrFstaDtFl,pthFldrFstaMirnsFls,pthFldrUTRFrmtdFlsSVMicro, \
	pthOutFldrSVMicro,pthToPhstCns,pthSVMicro,pthTmpFldr,pthToMAF, \
	pthToUCSCKEGG,sppRef,fOutDtSql,sqlFlSVMicro,lAttrbts)
	wrapper.addPrdctnslObjcts(lObjcts_trgt,lObjcts_toAdd, \
	['aMirnaSVMicroCnts','aMirnaSVMicroScrs'])
	#----------------------------
	#Ouput list of objects
	return lObjcts_trgt


########################################################
#~ Construct objects for sequences of 3'UTRs genes/isoforms, and lncRNAs 
# genes/isoforms from sql databases with annotations from all prediction 
# programs for humans.
########################################################
def mkl3UTRLncrnsRnAllPrdctns(pthFldrFstaUtr3cdgFl,fstaUtr3cdgFl, \
	rsltsUtr3cdgFldr,pthFldrFstaLncrnFl,fstaLncrnFl,rsltsLncrnFldr, \
	pthFldrFstaMirnsFls,fstasMirnaBsFl):
	"""
	Input: pthFldrFstaUtr3cdgFl is the path to the folder with the input 
	fasta file (with ENSEMBL header) for the 3'UTR of protein coding 
	genes/isoforms of interest. fstaUtr3cdgFl is the input fasta file 
	(with ENSEMBL header) for the 3'UTR of protein coding genes/isoforms 
	of interest. rsltsUtr3cdgFldr is the folder to write the results 
	for the 3'UTR of protein coding genes/isoforms of interest.
	pthFldrFstaLncrnFl is the path to the folder with the input fasta 
	file (with ENSEMBL header) for the lncRNA genes/isoforms of 
	interest. fstaLncrnFl is the input fasta file (with ENSEMBL header) 
	for the lncRNA genes/isoforms of interest. rsltsLncrnFldr is the 
	folder to write the results for the lncRNA genes/isoforms of 
	interest. pthFldrFstaMirnsFls is a path to the folder with the 
	sequence files of the objects. fstasMirnaBsFl is an array of miRNA 
	fasta-formated file names.
	Output: lObjctsUtr3cdg is a list of object with the results from all 
	the MRE predictions algorithms including in specific pointers for 
	3'UTR of protein coding genes/isoforms of interest included in 
	fstaUtr3cdgFl. lObjctsLncrn is a list of object with the resultsfrom 
	all the MRE predictions algorithms including in specific pointersfor 
	lncRNA genes/isoforms of interest included in fstaLncrnFl. 
	"""
	lObjctsUtr3cdg = mklObjctsRnAllPrdctns(pthFldrFstaUtr3cdgFl, \
	fstaUtr3cdgFl,rsltsUtr3cdgFldr,pthFldrFstaMirnsFls,fstasMirnaBsFl)
	lObjctsLncrn = mklObjctsRnAllPrdctns(pthFldrFstaLncrnFl, \
	fstaLncrnFl,rsltsLncrnFldr,pthFldrFstaMirnsFls,fstasMirnaBsFl)
	return lObjctsUtr3cdg,lObjctsLncrn


########################################################
#~ Return list of objects for target sequences of 3'UTRs genes/isoforms, 
# and lncRNAs genes/isoforms from background lists of objects.
########################################################
def rtrnTrgtlObjctsFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including in 
	specific pointers for 3'UTR of protein coding genes/isoforms of 
	interest. sSubstUtr3cdg is the subset of target 3'UTR sequence names 
	that must be present in lBckgrndUtr3cdg. lBckgrndLncrn is a list of 
	background object with the results from all the MRE predictions 
	algorithms including in specific pointers for lncRNA genes/isoforms 
	of interest. sSubstLncrn is the subset of target lncRNA sequence 
	names that must be present in lBckgrndLncrn.
	Output: lTrgtUtr3cdg is a copy of the object whose name are in 
	sSubstUtr3cdg for all pointers included in lBckgrndUtr3cdg list of 
	background objects. lTrgtLncrn is a copy of the object whose name 
	are in sSubstLncrn for all pointers included in lBckgrndLncrn listof 
	background objects.
	"""
	lTrgtUtr3cdg = __object__.mkArrySubStObjcts(lBckgrndUtr3cdg, \
	sSubstUtr3cdg)
	lTrgtLncrn = __object__.mkArrySubStObjcts(lBckgrndLncrn,sSubstLncrn)
	return lTrgtUtr3cdg,lTrgtLncrn


########################################################
#~ Run all the MRE prediction algorithms for a set miRNAs, background 
# and target sequences of 3'UTRs and lncRNAs.
########################################################
def rtrnLTrgtNLBckgrnd3UTRLncrns(pthFldrFstaUtr3cdgFl,fstaUtr3cdgFl, \
	rsltsUtr3cdgFldr,pthFldrFstaLncrnFl,fstaLncrnFl,rsltsLncrnFldr, \
	sSubstUtr3cdg,sSubstLncrn, \
	pthFldrFstaMirnsFls = '/exports/igmm/datastore/ponting-lab/CP008_BEDOYA-REINA_LNCRNA/GTeX_Mele.d/data.d/mirnas_data.d/hg_19.d', \
	fstasMirnaBsFl = 'mature_hg19.fa')
	"""
	Input: pthFldrFstaUtr3cdgFl is the path to the folder with the input 
	fasta file (with ENSEMBL header) for the 3'UTR of protein coding 
	genes/isoforms of interest. fstaUtr3cdgFl is the input fasta file 
	(with ENSEMBL header) for the 3'UTR of protein coding genes/isoforms 
	of interest. rsltsUtr3cdgFldr is the folder to write the results 
	for the 3'UTR of protein coding genes/isoforms of interest.
	pthFldrFstaLncrnFl is the path to the folder with the input fasta 
	file (with ENSEMBL header) for the lncRNA genes/isoforms of 
	interest. fstaLncrnFl is the input fasta file (with ENSEMBL header) 
	for the lncRNA genes/isoforms of interest. rsltsLncrnFldr is the 
	folder to write the results for the lncRNA genes/isoforms of 
	interest. sSubstUtr3cdg is the subset of target 3'UTR sequence names 
	that must be present in lBckgrndUtr3cdg. sSubstLncrn is the subset 
	of target lncRNA sequence names that must be present in fstaLncrnFl.
	pthFldrFstaMirnsFls is a path to the folder with the sequence files 
	of the objects. fstasMirnaBsFl is an array of miRNA fasta-formated 
	file names.
	Output: lBckgrndUtr3cdg is a list of object with the results fromall 
	the MRE predictions algorithms including in specific pointers for 
	3'UTR of protein coding genes/isoforms of interest included in 
	fstaUtr3cdgFl. lBckgrndLncrn is a list of object with the results
	from all the MRE predictions algorithms including in specific 
	pointers for lncRNA genes/isoforms of interest included in 
	fstaLncrnFl. lTrgtUtr3cdg is a copy of the object whose name are in 
	sSubstUtr3cdg for all pointers included in lBckgrndUtr3cdg list of 
	background objects. lTrgtLncrn is a copy of the object whose name 
	are in sSubstLncrn for all pointers included in lBckgrndLncrn listof 
	background objects.
	"""
	lBckgrndUtr3cdg,lBckgrndLncrn = mkl3UTRLncrnsRnAllPrdctns \
	(pthFldrFstaUtr3cdgFl,fstaUtr3cdgFl,rsltsUtr3cdgFldr, \
	pthFldrFstaLncrnFl,fstaLncrnFl,rsltsLncrnFldr,pthFldrFstaMirnsFls, \
	fstasMirnaBsFl)
	lTrgtUtr3cdg = __object__.mkArrySubStObjcts(lBckgrndUtr3cdg, \
	sSubstUtr3cdg)
	lTrgtLncrn = __object__.mkArrySubStObjcts(lBckgrndLncrn,sSubstLncrn)
	return lBckgrndUtr3cdg,lBckgrndLncrn,lTrgtUtr3cdg,lTrgtLncrn


#----------------------------
#Recipes to load or calculate correlation in expression and ECP
#----------------------------
########################################################
#~ Return the positions of genes/isoforms in an array of object from
# a subset of interest
########################################################
def rtrnPosSubstNm(lBckgrndObjcts,sSubstObjcts,vrbose=True):
	"""
	Input: lBckgrndObjcts is a list of background objects. sSubstObjcts 
	is a subset of names within lBckgrndObjcts. If vrbose is True the 
	log message is going to be printed.
	Output: aPosSubstNms is an array of the position for sSubstObjcts in 
	aNms.
	NOTE: only the interescting values are going to be returned.
	"""
	lTrgtObjcts = __object__.mkArrySubStObjcts(lBckgrndObjcts, \
	sSubstObjcts)
	dNmsPosInlBckgrndObjcts = dict([(objct.name,pos) for pos,objct in  \
	enumerate(lBckgrndObjcts)])
	aPosSubstNms = array([dNmsPosInlBckgrndObjcts[objct.name] for objct \
	in lTrgtObjcts])
	#----------------------------
	# report final results	
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('position for %s out of %s'%(len(aPosSubstNms),len(sSubstObjcts)),
		'nodes of interest were retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)	
	return aPosSubstNms,ovrAllogRun


########################################################
#~ Return a (sub)matrix with the values of two subsets from two list 
# of background objects (i.e. 3'UTR and lncRNAs) and a background matrix
# with values calculated between these lists of objects
########################################################
def rtrn3UTRLncrnsTrgAVlsFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,aVlsBckgrndUtr3cdgLncrn):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including in 
	specific pointers for 3'UTR of protein coding genes/isoforms of 
	interest. sSubstUtr3cdg is the subset of target 3'UTR sequence names 
	that must be present in lBckgrndUtr3cdg. lBckgrndLncrn is a list of 
	background object with the results from all the MRE predictions 
	algorithms including in specific pointers for lncRNA genes/isoforms 
	of interest. sSubstLncrn is the subset of target lncRNA sequence 
	names that must be present in lBckgrndLncrn. aVlsBckgrndUtr3cdgLncrn
	is a matrix with values between objects in lBckgrndUtr3cdg (in the
	same order), and objects in lBckgrndLncrn (in the same order).
	Output: aVlsTrgtUtr3cdgLncrn is a matrix for values for the names 
	of interest included in sSubstUtr3cdg and sSubstLncrn, taken from 
	the background matrix aVlsBckgrndUtr3cdgLncrn.
	NOTE: Values in aVlsBckgrndUtr3cdgLncrn are goin to be from the 
	sorted names in sSubstUtr3cdg and sSubstLncrn.
	"""
	assert len(lBckgrndUtr3cdg),len(lBckgrndLncrn) == \
	aVlsBckgrndUtr3cdgLncrn.shape
	aPosTrgtInBckgrndUtr3cdg = rtrnPosSubstNm(lBckgrndUtr3cdg, \
	sSubstUtr3cdg,vrbose=True)
	aPosTrgtInBckgrndLncrn = rtrnPosSubstNm(lBckgrndUtr3cdg, \
	sSubstUtr3cdg,vrbose=True)
	aVlsTrgtUtr3cdgLncrn = aVlsBckgrndUtr3cdgLncrn \
	[aPosTrgtInBckgrndUtr3cdg,:][:,aPosTrgtInBckgrndLncrn]
	lenaVlsTrgtUtr3cdg,lenaVlsTrgtLncrn = aVlsTrgtUtr3cdgLncrn.shape
	assert len(aPosTrgtInBckgrndUtr3cdg)==len(lenaVlsTrgtUtr3cdg) and \
	len(aPosTrgtInBckgrndLncrn)==len(lenaVlsTrgtLncrn) 
	return aVlsTrgtUtr3cdgLncrn
	

########################################################
#~ Return the positions for target sequences of 3'UTRs genes/isoforms, 
# and lncRNAs genes/isoforms from lists of objects and well as a 
# submatrix of values of interest.
########################################################
def rtrn3UTRLncrnsTrgtFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,aVlsBckgrndUtr3cdgLncrn):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	sSubstUtr3cdg is the subset of target 3'UTR sequence names that must 
	be present in lBckgrndUtr3cdg. lBckgrndLncrn is a list of background 
	object with the results from all the MRE predictions algorithms 
	including in specific pointers for lncRNA genes/isoforms ofinterest. 
	sSubstLncrn is the subset of target lncRNA sequence names that must 
	be present in lBckgrndLncrn. aVlsBckgrndUtr3cdgLncrn is a matrix 
	with values between objects in lBckgrndUtr3cdg (in the same order), 
	and objects in lBckgrndLncrn (in the same order).
	Output: lTrgtUtr3cdg is a copy of the object whose name are in 
	sSubstUtr3cdg for all pointers included in lBckgrndUtr3cdg list of 
	background objects. lTrgtLncrn is a copy of the object whose name 
	are in sSubstLncrn for all pointers included in lBckgrndLncrn listof 
	background objects. aVlsTrgtUtr3cdgLncrn is a matrix for values for 
	the names of interest included in sSubstUtr3cdg and sSubstLncrn, 
	taken from the background matrix aVlsBckgrndUtr3cdgLncrn.
	NOTE: Values in aVlsBckgrndUtr3cdgLncrn are goin to be from the 
	sorted names in sSubstUtr3cdg and sSubstLncrn.
	"""
	lTrgtUtr3cdg,lTrgtLncrn = rtrnTrgtlObjctsFrmBckgrnd(lBckgrndUtr3cdg, \
	sSubstUtr3cdg,lBckgrndLncrn,sSubstLncrn)
	aVlsTrgtUtr3cdgLncrn = rtrn3UTRLncrnsTrgAVlsFrmBckgrnd \
	(lBckgrndUtr3cdg,sSubstUtr3cdg,lBckgrndLncrn,sSubstLncrn, \
	aVlsBckgrndUtr3cdgLncrn)
	return lTrgtUtr3cdg,lTrgtLncrn,aVlsTrgtUtr3cdgLncrn
	

########################################################
#~ Load correlation in expression from a tsv file and the counting
# results from different predicting algorithms.
########################################################
def rnTSVCrrltnFlECPLncrnUtr3cdg(fCrrltnBckgrndLncrnUtr3cdgTsv, \
	lBckgrndUtr3cdg,lBckgrndLncrn,fldrOutECPPrws,rtrn3Mtrx=False):
	"""
	Input: fCrrltnBckgrndLncrnUtr3cdgTsv is a file in "tsv" format with 
	header have name of lncRNAs. In the same way, column 0 shall have 
	gene names. lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	lBckgrndLncrn is a list of background object with the results from 
	all the MRE predictions algorithms including in specific pointers 
	for lncRNA genes/isoforms ofinterest. fldrOutECPPrws is a temporal
	folder to write pairwise ecp values. Optionally, if rtrn3Mtrx is 
	True three matrices are going to be returned -instead of 1- 
	{aCrrltnBckgrndUtr3cdgLncrn for correlation in expression, 
	aECPBckgrndUtr3cdgLncrn for ecp values, and aPrctECPCrrltnUtr3cdgLncrn 
	for their product}.
	Output: aPrctECPCrrltnUtr3cdgLncrn is an array with the product of the
	average ECP value (determine by different programs) and the 
	Spearman's correlation in the values for the combinations of the 
	names in fCrrltnBckgrndLncrnUtr3cdgTsv. Optionally, if rtrn3Mtrx is 
	True three matrices are going to be returned -instead of 1- 
	{aCrrltnBckgrndUtr3cdgLncrn for correlation in expression, 
	aECPBckgrndUtr3cdgLncrn for ecp values, and aPrctECPCrrltnUtr3cdgLncrn 
	for their product}.
	"""	
	ovrAllogRun = []
	lPntrs = ['aMirnaRNAhCnts','aMirnaMirndCnts','aMirnaSVMicroCnts', \
	'aMirnaSVMicroCnts','aMirnaTrgtMnrCnts','aMirnaPITACnts', \
	'aMirnaMirMapCnts']
	#----------------------------
	#Load correlation in expression
	aCrrltnBckgrndUtr3cdgLncrn,aANmsBckgrndUtr3cdg,aANmsBckgrndLncrn = \
	coexpression.mkCrrltnFrmTSVfl(fCrrltnBckgrndLncrnUtr3cdgTsv)
	#----------------------------
	#Load ECP results	
	aMrnVlsBckgrndUtr3cdg = getattr(lBckgrndUtr3cdg,'aMirnaCnts')
	aMrnVlsBckgrndLncrn = getattr(lBckgrndLncrn,'aMirnaCnts')
	aECPBckgrndUtr3cdgLncrn = ecp.cmpECP(aMrnVlsBckgrndUtr3cdg, \
	aMrnVlsBckgrndLncrn,aANmsBckgrndUtr3cdg,aANmsBckgrndLncrn, \
	fldrOutECPPrws)
	for prdctPntr in lPntrs:
		aMrnVlsBckgrndUtr3cdg_toAdd = getattr(lBckgrndUtr3cdg,prdctPntr)
		aMrnVlsBckgrndLncrn_toAdd = getattr(lBckgrndLncrn,prdctPntr)
		aECPBckgrndUtr3cdgLncrn_toAdd = ecp.cmpECP \
		(aMrnVlsBckgrndUtr3cdg_toAdd,aMrnVlsBckgrndLncrn_toAdd, \
		aANmsBckgrndUtr3cdg,aANmsBckgrndLncrn,fldrOutECPPrws)	
		aECPBckgrndUtr3cdgLncrn,ovrAllogRun_toAdd = matrices. \
		rtrnIntrsctnTwoMtrx(aECPBckgrndUtr3cdgLncrn, \
		aECPBckgrndUtr3cdgLncrn_toAdd)
		ovrAllogRun.extend(ovrAllogRun_toAdd)
	aECPBckgrndUtr3cdgLncrn/=float32(len(lPntrs))
	#----------------------------
	#Obtain a product of correlation in expression and ECP
	aPrctECPCrrltnUtr3cdgLncrn,ovrAllogRun_toAdd = matrices. \
	rtrnIntrsctnTwoMtrx(aCrrltnBckgrndUtr3cdgLncrn, \
	aECPBckgrndUtr3cdgLncrn)
	ovrAllogRun.extend(ovrAllogRun_toAdd)
	if rtrn3Mtrx:
		return aCrrltnBckgrndUtr3cdgLncrn,aECPBckgrndUtr3cdgLncrn, \
		aPrctECPCrrltnUtr3cdgLncrn
	else:
		return aPrctECPCrrltnUtr3cdgLncrn


########################################################
#~ Load correlation in expression from a tsv file and the counting
# results from different predicting algorithms.
########################################################
def rnSqlCrrltnFlECPLncrnUtr3cdg(lBckgrndUtr3cdg,lBckgrndLncrn, \
	fldrOutECPPrws,lPntrs=['aMirnaRNAhCnts','aMirnaMirndCnts', \
	'aMirnaSVMicroCnts','aMirnaSVMicroCnts','aMirnaTrgtMnrCnts', \
	'aMirnaPITACnts','aMirnaMirMapCnts'],lSqlFls=[],rtrn3Mtrx=False):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	lBckgrndLncrn is a list of background object with the results from 
	all the MRE predictions algorithms including in specific pointers 
	for lncRNA genes/isoforms ofinterest. fldrOutECPPrws is a temporal
	folder to write pairwise ecp values. lPntrs is a list of pointers of
	interest to calculate ECP values. lSqlFls is a list of sql files 
	with the correlation in expression for tissues/samples of interest.
	Optionally, if rtrn3Mtrx is True three matrices are going to be 
	returned -instead of 1- {aCrrltnBckgrndUtr3cdgLncrn for correlation 
	in expression, aECPBckgrndUtr3cdgLncrn for ecp values, and 
	aPrctECPCrrltnUtr3cdgLncrn for their product}.
	Output: aPrctECPCrrltnUtr3cdgLncrn is an array with the product of the
	average ECP value (determine by different programs) and the 
	Spearman's correlation in the values for the combinations of the 
	names in lBckgrndUtr3cdg and lBckgrndLncrn. Optionally, if rtrn3Mtrx 
	is True three matrices are going to be returned -instead of 1- 
	{aCrrltnBckgrndUtr3cdgLncrn for correlation in expression, 
	aECPBckgrndUtr3cdgLncrn for ecp values, and aPrctECPCrrltnUtr3cdgLncrn 
	for their product}.
	"""	
	ovrAllogRun = []
	srtdUtr3cdgNms = [objct.name for objct in lBckgrndUtr3cdg]
	srtdLncrnNms = [objct.name for objct in lBckgrndLncrn]
	aCrrltnBckgrndUtr3cdgLncrn = coexpression.rtrnSqlFlCrrltn(lSqlFls[0], \
	srtdUtr3cdgNms,srtdLncrnNms,rtrnSgndCrrltn=True)#only positive
	#----------------------------
	#Load correlation in expression	results
	for sqlFl in lSqlFls[1:]:
		aCrrltnBckgrndUtr3cdgLncrn_toAdd = coexpression.rtrnSqlFlCrrltn \
		(sqlFl,srtdUtr3cdgNms,srtdLncrnNms,rtrnSgndCrrltn=True)#only positive
		aCrrltnBckgrndUtr3cdgLncrn,ovrAllogRun_toAdd = matrices. \
		rtrnIntrsctnTwoMtrx(aCrrltnBckgrndUtr3cdgLncrn, \
		aCrrltnBckgrndUtr3cdgLncrn_toAdd)
	aCrrltnBckgrndUtr3cdgLncrn/=float32(len(lSqlFls))#average coexpression
	#----------------------------
	#Load ECP results	
	aMrnVlsBckgrndUtr3cdg = getattr(lBckgrndUtr3cdg,'aMirnaCnts')
	aMrnVlsBckgrndLncrn = getattr(lBckgrndLncrn,'aMirnaCnts')
	aECPBckgrndUtr3cdgLncrn = ecp.cmpECP(aMrnVlsBckgrndUtr3cdg, \
	aMrnVlsBckgrndLncrn,aANmsBckgrndUtr3cdg,aANmsBckgrndLncrn, \
	fldrOutECPPrws)
	for prdctPntr in lPntrs:
		aMrnVlsBckgrndUtr3cdg_toAdd = getattr(lBckgrndUtr3cdg,prdctPntr)
		aMrnVlsBckgrndLncrn_toAdd = getattr(lBckgrndLncrn,prdctPntr)
		aECPBckgrndUtr3cdgLncrn_toAdd = ecp.cmpECP \
		(aMrnVlsBckgrndUtr3cdg_toAdd,aMrnVlsBckgrndLncrn_toAdd, \
		aANmsBckgrndUtr3cdg,aANmsBckgrndLncrn,fldrOutECPPrws)	
		aECPBckgrndUtr3cdgLncrn,ovrAllogRun_toAdd = matrices. \
		rtrnIntrsctnTwoMtrx(aECPBckgrndUtr3cdgLncrn, \
		aECPBckgrndUtr3cdgLncrn_toAdd)
		ovrAllogRun.extend(ovrAllogRun_toAdd)
	aECPBckgrndUtr3cdgLncrn/=float32(len(lPntrs))
	#----------------------------
	#Obtain a product of correlation in expression and ECP
	aPrctECPCrrltnUtr3cdgLncrn,ovrAllogRun_toAdd = matrices. \
	rtrnIntrsctnTwoMtrx(aCrrltnBckgrndUtr3cdgLncrn, \
	aECPBckgrndUtr3cdgLncrn)
	ovrAllogRun.extend(ovrAllogRun_toAdd)
	if rtrn3Mtrx:
		return aCrrltnBckgrndUtr3cdgLncrn,aECPBckgrndUtr3cdgLncrn, \
		aPrctECPCrrltnUtr3cdgLncrn
	else:
		return aPrctECPCrrltnUtr3cdgLncrn


########################################################
#~ Obtain target data information from a background dataset
########################################################
def rtrnTrgtBckrngVrblsFrmSql(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,fldrOutECPPrws,
	lPntrs=['aMirnaRNAhCnts','aMirnaMirndCnts', \
	'aMirnaSVMicroCnts','aMirnaSVMicroCnts','aMirnaTrgtMnrCnts', \
	'aMirnaPITACnts','aMirnaMirMapCnts'],lSqlFls=[],rtrn3Mtrx=True):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	sSubstUtr3cdg is the subset of target 3'UTR sequence names that must 
	be present in lBckgrndUtr3cdg. lBckgrndLncrn is a list of background 
	object with the results from all the MRE predictions algorithms 
	including in specific pointers for lncRNA genes/isoforms ofinterest. 
	sSubstLncrn is the subset of target lncRNA sequence names that must 
	be present in lBckgrndLncrn. fldrOutECPPrws is a temporal folder to 
	write pairwise ecp values. lPntrs is a list of pointers of interest 
	to calculate ECP values. lSqlFls is a list of sql files with the 
	correlation in expression for tissues/samples of interest.
	Optionally, if rtrn3Mtrx is True three matrices are going to be 
	returned -instead of 1- {aCrrltnBckgrndUtr3cdgLncrn for correlation 
	in expression, aECPBckgrndUtr3cdgLncrn for ecp values, and 
	aPrctECPCrrltnVlsAVlsB for their product}.
	Output: lTrgtUtr3cdg is a copy of the object whose name are in 
	sSubstUtr3cdg for all pointers included in lBckgrndUtr3cdg list of 
	background objects. lTrgtLncrn is a copy of the object whose name 
	are in sSubstLncrn for all pointers included in lBckgrndLncrn listof 
	background objects. aPrctECPCrrltnBckgrndUtr3cdgLncrn is an array 
	with the product of the average ECP value (determine by different 
	programs) and the Spearman's correlation in the values for the 
	combinations of the names in lBckgrndUtr3cdg and lBckgrndLncrn. 
	aCrrltnBckgrndUtr3cdgLncrn is the same matrix for correlation in 
	expression, and aECPBckgrndUtr3cdgLncrn for ecp values. 
	aCrrltnTrgtUtr3cdgLncrn, aECPTrgtUtr3cdgLncrn and 
	aPrctECPCrrltnTrgtUtr3cdgLncrn are the same for the target objects.
	"""
	aCrrltnBckgrndUtr3cdgLncrn,aECPBckgrndUtr3cdgLncrn, \
	aPrctECPCrrltnBckgrndUtr3cdgLncrn = rnSqlCrrltnFlECPLncrnUtr3cdg \
	(lBckgrndUtr3cdg,lBckgrndLncrn,fldrOutECPPrws,lPntrs,lSqlFls, \
	rtrn3Mtrx)
	#----------------------------
	#Obtain values for target correlations
	lTrgtUtr3cdg,lTrgtLncrn,aCrrltnTrgtUtr3cdgLncrn = \
	rtrn3UTRLncrnsTrgtFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,aCrrltnBckgrndUtr3cdgLncrn)
	#----------------------------
	#Obtain values for target ECP values
	lTrgtUtr3cdgA,lTrgtLncrnA,aECPTrgtUtr3cdgLncrn = \
	rtrn3UTRLncrnsTrgtFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,aECPBckgrndUtr3cdgLncrn)
	assert lTrgtUtr3cdg==lTrgtUtr3cdgA and lTrgtLncrn==lTrgtLncrnA
	#----------------------------
	#Obtain values for target correlations*ECP
	lTrgtUtr3cdgA,lTrgtLncrnA,aPrctECPCrrltnTrgtUtr3cdgLncrn = \
	rtrn3UTRLncrnsTrgtFrmBckgrnd(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,aPrctECPCrrltnBckgrndUtr3cdgLncrn)
	assert lTrgtUtr3cdg==lTrgtUtr3cdgA and lTrgtLncrn==lTrgtLncrnA
	#----------------------------
	return aCrrltnBckgrndUtr3cdgLncrn,aECPBckgrndUtr3cdgLncrn, \
	aPrctECPCrrltnBckgrndUtr3cdgLncrn,lTrgtUtr3cdg,lTrgtLncrn, \
	aCrrltnTrgtUtr3cdgLncrn,aECPTrgtUtr3cdgLncrn, \
	aPrctECPCrrltnTrgtUtr3cdgLncrn


########################################################
#~ Return an statistic over an array of array of statiscs on miRNA of 
# 3'UTRs genes/isoforms and lncRNAs genes/isoforms.
########################################################
def rtrnStstMirnArrys(lObjctUtr3cdg,lObjctLncrn,lPntrs= \
	['aMirnaRNAhCnts','aMirnaMirndCnts','aMirnaSVMicroCnts', \
	'aMirnaSVMicroCnts','aMirnaTrgtMnrCnts','aMirnaPITACnts', \
	'aMirnaMirMapCnts']):
	"""
	Input: lObjctUtr3cdg is the list of objects with array of miRNA 
	metrics of interest in lPntrs for 3'UTR objects. lObjctLncrn is the 
	list of objects with array of miRNA metrics of interest in lPntrs 
	for lncRNAs. lPntrs is the list of pointer with values for MREs
	predicted from different algorithms.
	Output: aMirnMetrcMskd is an array with array with the average 
	MREs values caluclated for different algorithms.
	"""
	aUtr3cdgMirnMetrc = array([getattr(dt,lPntrs[0]) for dt in \
	lObjctUtr3cdg])
	aLncrnMirnMetrc = array([getattr(dt,lPntrs[0]) for dt in \
	lObjctLncrn])
	aMirnMetrcMskd = randomz.cmptMirnaStatRowClmnPrs(aUtr3cdgMirnMetrc, \
	aLncrnMirnMetrc)
	for pntr in lPntrs[1:]:
		aUtr3cdgMirnMetrc = array([getattr(dt,pntr) for dt in \
		lObjctUtr3cdg])
		aLncrnMirnMetrc = array([getattr(dt,pntr) for dt in lObjctLncrn])
		aMirnMetrcMskd_toAdd = randomz.cmptMirnaStatRowClmnPrs \
		(aUtr3cdgMirnMetrc,aLncrnMirnMetrc)
		aMirnMetrcMskd+=aMirnMetrcMskd_toAdd
	aMirnMetrcMskd/=float32(len(lPntrs))
	return aMirnMetrcMskd
	

#----------------------------
#Recipes to test competiton
#----------------------------

########################################################
#~ Calculate model of length distribution for two list of target objects 
# (i.e. 3'UTR and lncRNAs)
########################################################
def rtrnLenMdlUtr3cdgNLncrn(lTrgtUtr3cdg,outLenMdlPltUtr3cdgTrnscrpts, \
	lTrgtLncrn,outLenMdlPltLncrnaTrnscrpts):
	"""
	Input: lTrgtUtr3cdg is a list of target object for 3'UTRs of protein
	coding genes. outLenMdlPltUtr3cdgTrnscrpts is the complete path to 
	a file to save the model plot for lTrgtUtr3cdg. lTrgtLncrn is a list 
	of target object for lncRNAs. outLenMdlPltLncrnaTrnscrpts is the 
	complete path to a file to save the model plot for lTrgtLncrn.
	Output: ftdMdlPrmtrsLenONLYUtr3cdg is the model for length 
	distribution of objects in lTrgtUtr3cdg. ftdMdlPrmtrsLenLncrna is 
	the model for length distribution of objects in lTrgtLncrn.
	"""
	ftdMdlPrmtrsLenLncrna,messgMdlLenLncrna = \
	sampling.fitMdl(lTrgtLncrn,outLenMdlPltLncrnaTrnscrpts)
	ftdMdlPrmtrsLenONLYUtr3cdg,messgLenONLYUtr3cdg = sampling. \
	fitMdl(lTrgtUtr3cdg,outLenMdlPltUtr3cdgTrnscrpts)
	return ftdMdlPrmtrsLenONLYUtr3cdg,ftdMdlPrmtrsLenLncrna


########################################################
#~ Calculate significance of ECP for two list of target objects (i.e. 
# 3'UTR and lncRNAs) with positive correlation in expression
########################################################
def rtrnSgnfncUtr3cdgNLncrnCmpttn(lBckgrndUtr3cdg,sSubstUtr3cdg, \
	lBckgrndLncrn,sSubstLncrn,fldrOutECPPrws,lPntrs=['aMirnaRNAhCnts', \
	'aMirnaMirndCnts','aMirnaSVMicroCnts','aMirnaSVMicroCnts', \
	'aMirnaTrgtMnrCnts','aMirnaPITACnts','aMirnaMirMapCnts'],lSqlFls=[], \
	rtrn3Mtrx=True, \
	outLenMdlPltUtr3cdgTrnscrpts, \#output file for plot on 3'UTR lengths
	outLenMdlPltLncrnaTrnscrpts, \#output file for plot on lncRNA lengths
	outRndmztnFldr, \#output folder for randomizations
	outMirnStstFl, \#output folder for mirna stats file
	outPltFl):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	sSubstUtr3cdg is the subset of target 3'UTR sequence names that must 
	be present in lBckgrndUtr3cdg. lBckgrndLncrn is a list of background 
	object with the results from all the MRE predictions algorithms 
	including in specific pointers for lncRNA genes/isoforms ofinterest. 
	sSubstLncrn is the subset of target lncRNA sequence names that must 
	be present in lBckgrndLncrn. fldrOutECPPrws is a temporal folder to 
	write pairwise ecp values. lPntrs is a list of pointers of interest 
	to calculate ECP values. lSqlFls is a list of sql files with the 
	correlation in expression for tissues/samples of interest.
	Optionally, if rtrn3Mtrx is True three matrices are going to be 
	returned -instead of 1- {aCrrltnBckgrndUtr3cdgLncrn for correlation 
	in expression, aECPBckgrndUtr3cdgLncrn for ecp values, and 
	aPrctECPCrrltnVlsAVlsB for their product}. 
	outLenMdlPltUtr3cdgTrnscrpts is the output file to write the plot 
	for the distribution and model of the length of 3'UTRs. 
	outLenMdlPltLncrnaTrnscrpts is the output file to write the plot for 
	the distribution and model of the length of lncRNAs. outRndmztnFldr
	is the folder to write the randomization results. outMirnStstFl is 
	the file to write the enrichement of miRNAs. outPltFl is the file to
	write the output plot for the significance of the target metric in 
	the background randomizations.
	Output: ovrAllogRun is the log run for a full randomization. p_val 
	is the p-value for the significance of the target metric given a 
	full randomization of data. ovrAllogRunUtr3cdg is the log run for 
	the randomization of Utr3cdg. p_valUtr3cdg is the p-value for the 
	significance of the target metric given a full randomization of 
	Utr3cdg. ovrAllogRunLncrn is the log run for the randomization of 
	lncRNAs. p_valLncrn is the p-value for the significance of the 
	target metric given a full randomization of lncRNAs.
	"""
	#----------------------------
	#Set start variables
	aMirnNms = getattr(lBckgrndUtr3cdg[0],'aMirnaNms')
	mskRowClmnDflt=None
	statstc='mean'
	aPosClmnsToGns=None
	aPosRowsToGns=None
	mirnDtype='cnt'
	#----------------------------
	#return background and target variables
	aCrrltnBckgrndUtr3cdgLncrn,aECPBckgrndUtr3cdgLncrn, \
	aPrctECPCrrltnBckgrndUtr3cdgLncrn,lTrgtUtr3cdg,lTrgtLncrn, \
	aCrrltnTrgtUtr3cdgLncrn,aECPTrgtUtr3cdgLncrn, \
	aPrctECPCrrltnTrgtUtr3cdgLncrn =  rtrnTrgtBckrngVrblsFrmSql \
	(lBckgrndUtr3cdg,sSubstUtr3cdg,lBckgrndLncrn,sSubstLncrn, \
	fldrOutECPPrws,lPntrs,lSqlFls,rtrn3Mtrx)
	#----------------------------
	#model length of target data
	ftdMdlPrmtrsLenONLYUtr3cdg,ftdMdlPrmtrsLenLncrna = \
	rtrnLenMdlUtr3cdgNLncrn(lTrgtUtr3cdg,outLenMdlPltUtr3cdgTrnscrpts, \
	lTrgtLncrn,outLenMdlPltLncrnaTrnscrpts)
	#----------------------------
	#obtain target dataset
	mskNgtvCrrltnsTrgtUtr3cdgLncrn = ma.masked_invalid \
	(aCrrltnTrgtUtr3cdgLncrn).mask
	aTrgtMetrcMskd = ma.array(aECPTrgtUtr3cdgLncrn,mask = \
	mskNgtvCrrltnsTrgtUtr3cdgLncrn.mask)
	aTrgtMetrcMskd.fill_value = nan
	aTrgtMetrcMskd = ma.masked_invalid(aTrgtMetrcMskd.filled())
	#----------------------------
	#obtain target dataset for miRNAs
	aTrgtMirnMetrcMskd = rtrnStstMirnArrys(lTrgtUtr3cdg,lTrgtLncrn, \
	lPntrs)
	aTrgtMirnMetrcMskd.fill_value = nan
	aTrgtMirnMetrcMsk = aTrgtMetrcMskd.mask[:,:,newaxis]*ones \
	(aTrgtMirnMetrcMskd.shape)#mask for target mirnas
	aTrgtMirnMetrcMskd.mask = ma.mask_or(aTrgtMirnMetrcMsk.mask, \
	aTrgtMirnMetrcMsk)
	#----------------------------
	#obtain background dataset
	mskNgtvCrrltnsBckgrndUtr3cdgLncrn = ma.masked_invalid \
	(aCrrltnBckgrndUtr3cdgLncrn).mask
	aBckgrndMetrcMskd = ma.array(aECPBckgrndUtr3cdgLncrn,mask = \
	mskNgtvCrrltnsBckgrndUtr3cdgLncrn.mask)
	aBckgrndMetrcMskd.fill_value = nan
	aBckgrndMetrcMskd = ma.masked_invalid(aBckgrndMetrcMskd.filled())	
	#----------------------------
	#obtain background dataset for miRNAs
	aBckgrndMirnMetrcMskd = rtrnStstMirnArrys(lBckgrndUtr3cdg, \
	lBckgrndLncrn,lPntrs)
	aBckgrndMirnMetrcMskd.fill_value = nan
	aBckgrndMirnMetrcMsk = aBckgrndMetrcMskd.mask[:,:,newaxis]*ones \
	(aBckgrndMirnMetrcMskd.shape)#mask for target mirnas
	aBckgrndMirnMetrcMskd.mask = ma.mask_or(aBckgrndMirnMetrcMsk.mask, \
	aBckgrndMirnMetrcMsk)
	#----------------------------
	#run full randomizations
	lBckgrndUtr3cdgLens = [dt.len for dt in lBckgrndUtr3cdg]
	lBckgrndLncrnLens = [dt.len for dt in lBckgrndLncrn]
	ovrAllogRun,p_val = randomz.wrprFllRndmztn(aBckgrndMetrcMskd, \
	lBckgrndUtr3cdgLens,lBckgrndLncrnLens,ftdMdlPrmtrsLenONLYUtr3cdg, \
	ftdMdlPrmtrsLenLncrna,aTrgtMetrcMskd=aTrgtMetrcMskd,outPltFl= \
	outPltFl,statstc=statstc,aPosClmnsToGns=aPosClmnsToGns, \
	aPosRowsToGns=aPosRowsToGns,aBckgrndMirnMetrcMskd= \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd=aTrgtMirnMetrcMskd, \
	mirnDtype=mirnDtype,outMirnStstFl=outMirnStstFl,aMirnNms=aMirnNms)	
	#----------------------------
	#run 3'UTR randomizations
	outPltFl,extsn = os.path.splitext()
	outPltFl = '.'.join(['%s_Utr3cdg'%outPltFl,extsn])
	ovrAllogRunUtr3cdg,p_valUtr3cdg = \
	wrprRowRndmztn(aBckgrndMetrcMskd,lBckgrndUtr3cdgLens, \
	ftdMdlPrmtrsLenLncrna,aTrgtMetrcMskd=aTrgtMetrcMskd, \
	outPltFl=outPltFl,statstc=statstc,aMirnNms=aMirnNms)
	#----------------------------
	#run lncRNA randomizations
	outPltFl,extsn = os.path.splitext()
	outPltFl = '.'.join(['%s_Lncrn'%outPltFl,extsn])
	ovrAllogRunLncrn,p_valLncrn = \
	wrprClmRndmztn(aBckgrndMetrcMskd,lBckgrndLncrnLens, \
	ftdMdlPrmtrsLenONLYUtr3cdg,aTrgtMetrcMskd=aTrgtMetrcMskd, \
	outPltFl=outPltFl,statstc=statstc,aMirnNms=aMirnNms)
	return ovrAllogRun,p_val,ovrAllogRunUtr3cdg,p_valUtr3cdg, \
	ovrAllogRunLncrn,p_valLncrn


#----------------------------
#Recipes for graph analysis
#----------------------------

########################################################
#~ Make a bi-partite graph out of two loaded matrices, intersect them, 
# assess the significance of the three matrices, and output gml files 
# for the three of them.
########################################################
def mkGrphTstSgnfcnc(aVlsAVlsB_1,aVlsAVlsB_2,aVlsANms,aVlsBNms, \
	outGMLFl_1,outGMLFl_2,outGMLFl_I,vrbse=True):
	"""
	Input: Two different matrices with the same dimensions {aVlsAVlsB1,
	aVlsAVlsB2}. aVlsANms is an array with names in the same order as 
	rows in aVlsAVlsB_1 and aVlsAVlsB_2. aVlsBIsfNms is an array with 
	names in the same order as columns in aVlsAVlsB_1 and aVlsAVlsB_2.  
	outGMLFl_1  is the output file in GML format with a subgraph of 
	nodes annotated with red color (those edges interrsecting with 
	outGMLFl_2). outGMLFl_2 is the output file in GML format with a 
	subgraph of nodes annotated with red color (those edges 
	interrsecting with outGMLFl_1). outGMLFl_I in the graph intersecting 
	aVlsAVlsB_1 and aVlsAVlsB_2 whose weight is the product of the 
	weight of their edges. If vrbse is True, the log for the run is 
	going to be printed.
	Output: ovrAllogRun is the log of the run.
	NOTE: Program makes outGMLFl_1, outGMLFl_2, and outGMLFl_I. 
	outGMLFl_1 is the output file in GML format with a subgraph of nodes 
	annotated with red color (those edges interrsecting without GMLFl_2). 
	outGMLFl_2 is the output file in GML format with a subgraph of nodes 
	annotated with red color (those edges interrsecting with outGMLFl_1). 
	outGMLFl_I in the graph intersecting aVlsAVlsB_1 and aVlsAVlsB_2 
	whose weight is the product of the weight of their edges. 
	"""
	#----------------------------
	#~ Make graphs from input matrices
	ovrAllogRun = []
	G_1,ovrAllogRun_1 = networks.mkNtwrkXFrmMtrx(aVlsAVlsB_1,aVlsANms, \
	aVlsBNms,vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_1)
	G_2,ovrAllogRun_2 = networks.mkNtwrkXFrmMtrx(aVlsAVlsB_2,aVlsANms, \
	aVlsBNms,vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_2)
	#----------------------------
	#~ Run randomizations
	ovrAllogRun_rndmztn = statistics.wrpRndmzNtwrk(G,n_rndmStps=1000, \
	vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_rndmztn)
	#----------------------------
	#~ Write GML files
	ovrAllogRun_I = anntNtwrkXEdgsNotInRef(G_1,G_2,vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_I)
	networks.writeNtwrkXToGML(G_1,outGMLFl_1)
	networks.writeNtwrkXToGML(G_2,outGMLFl_2)
	#----------------------------
	#~ Retrieve network intersection
	aVlsAVlsBI,ovrAllogRun_I = rtrnIntrsctnTwoMtrx(aVlsAVlsB1, \
	aVlsAVlsB2,vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_I)
	G_I,ovrAllogRun_I = networks.mkNtwrkXFrmMtrx(aVlsAVlsBI,aVlsANms, \
	aVlsBNms,vrbse=vrbse)
	ovrAllogRun.extend(ovrAllogRun_I)
	networks.writeNtwrkXToGML(G_I,outGMLFl_I)
	return ovrAllogRun


########################################################
#~ Build a network for a target matrix and assess it significance
########################################################
def mkNtrkFrmMtrxNSgnfcnc(lBckgrndUtr3cdg,lBckgrndLncrn,sSubstUtr3cdg, \
	sSubstLncrn,fldrOutECPPrws, \
	,lPntrs=['aMirnaRNAhCnts','aMirnaMirndCnts','aMirnaSVMicroCnts', \
	'aMirnaSVMicroCnts','aMirnaTrgtMnrCnts','aMirnaPITACnts', \
	'aMirnaMirMapCnts'], \
	lSqlFls=[], \
	G_targetGMLOut, \
	G_targetGMLAnntdOut):
	"""
	Input: lBckgrndUtr3cdg is a list of background object with the 
	results from all the MRE predictions algorithms including inspecific 
	pointers for 3'UTR of protein coding genes/isoforms of interest. 
	lBckgrndLncrn is a list of background object with the results from 
	all the MRE predictions algorithms including in specific pointers 
	for lncRNA genes/isoforms ofinterest. sSubstUtr3cdg is the subset of 
	target 3'UTR sequence names that must be present in lBckgrndUtr3cdg. 
	sSubstLncrn is the subset of target lncRNA sequence names that must 
	be present in lBckgrndLncrn. fldrOutECPPrws is a temporal folder to 
	write pairwise ecp values. lPntrs is a list of pointers of interest 
	to calculate ECP values. lSqlFls is a list of sql files with the 
	correlation in expression for tissues/samples of interest.
	G_targetGMLOut is the path to the GML output network. 
	G_targetGMLAnntdOut is the path to the annotated GML output network. 
	Output: G_targetGMLOut is the path to the GML output network with 
	empirical frequencies of each edge in the target graph. 
	G_targetGMLAnntdOut is the path to the annotated version (i.e. with
	common names) of G_targetGMLOut. ovrAllogRun is the log run.
	"""
	#----------------------------
	#~ Calculate background product of ECP and correlations for 
	# background
	ovrAllogRun = []
	aPrctECPCrrltnBckgrndUtr3cdgLncrn = rnSqlCrrltnFlECPLncrnUtr3cdg \
	(lBckgrndUtr3cdg,lBckgrndLncrn,fldrOutECPPrws,lPntrs,lSqlFls, \
	rtrn3Mtrx=False)
	#----------------------------
	#~ Retrieve target list of obecjts
	lTrgtUtr3cdg,lTrgtLncrn = rtrnTrgtlObjctsFrmBckgrnd(lBckgrndUtr3cdg, \
	sSubstUtr3cdg,lBckgrndLncrn,sSubstLncrn)
	#----------------------------
	#~ Retrieve the array for the target product of ECP and correlations
	aPrctECPCrrltnTrgtUtr3cdgLncrn = rtrn3UTRLncrnsTrgAVlsFrmBckgrnd \
	(lBckgrndUtr3cdg,sSubstUtr3cdg,lBckgrndLncrn,sSubstLncrn, \
	aPrctECPCrrltnBckgrndUtr3cdgLncrn)
	#----------------------------
	#~ Make graph
	aVlsUtr3cdgNms = array([dt.name for dt in lTrgtUtr3cdg])
	aVlsLncrnNms = array([dt.name  for dt in lTrgtLncrn])
	G,ovrAllogRun_ntwrk = networks.mkNtwrkXFrmMtrx \
	(aPrctECPCrrltnTrgtUtr3cdgLncrn,aVlsUtr3cdgNms,aVlsLncrnNms)
	ovrAllogRun.extend(ovrAllogRun_ntwrk)
	#----------------------------
	#~ Make randomizations on graph
	ovrAllogRun_sts = statistics.wrpRndmzNtwrk(G)
	ovrAllogRun.extend(ovrAllogRun_sts)
	#----------------------------
	#~ Make GML file
	writeNtwrkXToGML(G,G_targetGMLOut)
	#----------------------------
	#~ Annotate graph with common names
	dLbldNmAnntn = dict([(dt.name,dt.cmmNm) for dt in lTrgtUtr3cdg])
	dLbldNmAnntn.update(dict([(dt.name,dt.cmmNm) for dt in lTrgtLncrn]))
	anntGMLGrph(dLbldNmAnntn,G_targetGMLOut,G_targetGMLAnntdOut)
	return ovrAllogRun

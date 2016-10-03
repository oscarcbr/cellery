#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_objects.py part of cellery (ceRNAs linking inference)
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
Test the building objects-related methods
"""

########################################################
#~ Import libraries.
########################################################
from cellery import __db__
from cellery import __object__
import os


########################################################
#~ Download isoform databases -execute-.
########################################################
def dwnldIsfrmDB(sTrnscrpts,sTrnscrpts3UTRs,sLncrnasTrnscrpts, \
	fOutPCTrscrpts,fOutPCTrscrptsSbst,fOutUtr3cdgTrnscrpts, \
	fOutUtr3cdgTrnscrptsSbst,fOutLncrnaTrnscrpts):
	"""
	"""
	#~ Retrive from ENSEMBL database for all protein coding isoforms.
	if not os.path.exists(fOutPCTrscrpts):
		__db__.rtrnPCTrscrpts(fOutPCTrscrpts)
	#~ Retrive from ENSEMBL database from a set of protein coding 
	# isoforms.
	if not os.path.exists(fOutPCTrscrptsSbst):
		__db__.rtrnPCTrscrptsSubSt(sTrnscrpts,fOutPCTrscrptsSbst)
	#~ Retrive from ENSEMBL database 3'UTRs from all protein coding 
	# isoforms.
	if not os.path.exists(fOutUtr3cdgTrnscrpts):
		__db__.rtrn3UTRhg19Trnscrpts(fOutUtr3cdgTrnscrpts)
	#~ Retrive from ENSEMBL database 3'UTRs from a set of protein coding 
	# isoforms of interest.
	if not os.path.exists(fOutUtr3cdgTrnscrptsSbst):
		__db__.rtrn3UTRhg19TrnscrptsSubSt(sTrnscrpts3UTRs, \
		fOutUtr3cdgTrnscrptsSbst)
	#~ Retrive from ENSEMBL database lncRNAs from a set of protein  
	# coding isoforms of interest.
	if not os.path.exists(fOutLncrnaTrnscrpts):
		__db__.rtrnLncRNATrnscrptsSubSt(sLncrnasTrnscrpts, \
		fOutLncrnaTrnscrpts)
	return 0
	

########################################################
#~ Download gene databases -execute-.
########################################################
def dwnldGnDB(sGns,sGns3UTRs,sLncrnas,fOutPCGenes,fOutPCGenesSbst, \
	fOutUtr3cdg,fOutUtr3cdgSbst,fOutLncrna):
	"""
	"""
	#~ Retrive from ENSEMBL database for all protein coding genes.
	if not os.path.exists(fOutPCGenes):
		__db__.rtrnPCGenes(fOutPCGenes)
	#~ Retrive from ENSEMBL database from a set of protein coding genes.
	if not os.path.exists(fOutPCGenesSbst):
		__db__.rtrnPCGenesSubSt(sGns,fOutPCGenesSbst)
	#~ Retrive from ENSEMBL database 3'UTRs from all protein coding genes.
	if not os.path.exists(fOutUtr3cdg):
		__db__.rtrn3UTRGnsHg19(fOutUtr3cdg)
	#~ Retrive from ENSEMBL database 3'UTRs from a set of protein coding 
	# genes of interest.
	if not os.path.exists(fOutUtr3cdgSbst):
		__db__.rtrn3UTRGnsHg19SubSt(sGns3UTRs,fOutUtr3cdgSbst)
	#~ Retrive from ENSEMBL database lncRNAs from a set of protein coding 
	# genes of interest.
	if not os.path.exists(fOutLncrna):
		__db__.rtrnLncRNAGnsSubSt(sLncrnas,fOutLncrna)
	return 0


########################################################
#~ Build initial database for sequences of protein coding genes and 
# lncRNAs.
########################################################
#~ Database output folder
dbFldr = '/gpfs/igmmfs01/eddie/ponting-lab/cellery/db/GTeX_Human.d'
#--------------------------------------------------------
#~ Download GENCODE lncRNA isoform database (August 1, 2016 -release 25-)
if not os.path.exists(os.path.join(dbFldr, \
	'gencode.v19.lncRNA_transcripts.fa')):
	os.system('wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz -P %s' \
	%dbFldr)
	os.system('gunzip %s'%os.path.join(dbFldr, \
	'gencode.v19.lncRNA_transcripts.fa.gz'))
#~ Download GENCODE PCT isoform database (August 1, 2016 -release 25-)
if not os.path.exists(os.path.join(dbFldr, \
	'gencode.v19.pc_transcripts.fa')):
	os.system('wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz -P %s' \
	%dbFldr)
	os.system('gunzip %s'%os.path.join(dbFldr, \
	'gencode.v19.pc_transcripts.fa.gz'))
#--------------------------------------------------------
#~ All protein coding isoforms -file-
fOutPCTrscrpts = os.path.join(dbFldr,'hg19_ENSEMBL_AllPCT_Isfrms.fas')
#~ Subset of protein coding isoforms -All annotated by GENCODE- -file-
fOutPCTrscrptsSbst = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODEPCT_Isfrms.fas')
#~ All 3'UTRs from protein coding isoforms -file-
fOutUtr3cdgTrnscrpts = os.path.join(dbFldr, \
'hg19_ENSEMBL_All3UTR_Isfrms.fas')
#~ Subset of 3'UTRs from protein coding isoforms -file-
fOutUtr3cdgTrnscrptsSbst = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODE3UTR_Isfrms.fas')
#~ Subset of lncRNA isoforms -All annotated by GENCODE- -file-
fOutLncrnaTrnscrpts = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODELncrns_Isfrms.fas')
#~ All protein coding genes -file-
fOutPCGenes = os.path.join(dbFldr,'hg19_ENSEMBL_AllPCT_Genes.fas')
#~ Subset of protein coding genes -file-
fOutPCGenesSbst = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODEPCT_Genes.fas')
#~ All 3'UTRs from protein coding genes -file-
fOutUtr3cdg = os.path.join(dbFldr, \
'hg19_ENSEMBL_All3UTR_Genes.fas')
#~ Subset of 3'UTRs from protein coding genes -file-
fOutUtr3cdgSbst = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODE3UTR_Genes.fas')
#~ Subset of lncRNA genes -All annotated by GENCODE- -file-
fOutLncrna = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODELncrns_Genes.fas')
#--------------------------------------------------------
#~ Subset of protein coding isoforms -All annotated by GENCODE-
fSbstPCT = os.path.join(dbFldr,'gencode.v19.pc_transcripts.fa')
sTrnscrpts = set([l.split('.')[0][1:] for l in open(fSbstPCT) if \
l[0]=='>'])
#~ Subset of 3'UTRs from protein coding isoforms -All annotated by 
# GENCODE-
sTrnscrpts3UTRs = set([l.split('.')[0][1:] for l in open(fSbstPCT) if \
l[0]=='>'])
#~ Subset of lncRNA isoforms -All annotated by GENCODE-
fSbstLncrns = os.path.join(dbFldr, \
'gencode.v19.lncRNA_transcripts.fa')
sLncrnasTrnscrpts = set([l.split('.')[0][1:] for l in open(fSbstLncrns) \
if l[0]=='>'])
#~ Subset of protein coding genes -All annotated by GENCODE-
sGns = set([l.split('|')[1].split('.')[0] for l in open(fSbstPCT) if \
l[0]=='>'])
#~ Subset of 3'UTRs from protein coding genes -All annotated by GENCODE-
sGns3UTRs = set([l.split('|')[1].split('.')[0] for l in open(fSbstPCT) \
if l[0]=='>'])
#~ Subset of lncRNA genes -All annotated by GENCODE-
sLncrnas = set([l.split('|')[1].split('.')[0] for l in open(fSbstLncrns) \
if l[0]=='>'])
#--------------------------------------------------------
#~ Execute
dwnldIsfrmDB(sTrnscrpts,sTrnscrpts3UTRs,sLncrnasTrnscrpts, \
fOutPCTrscrpts,fOutPCTrscrptsSbst,fOutUtr3cdgTrnscrpts, \
fOutUtr3cdgTrnscrptsSbst,fOutLncrnaTrnscrpts)
dwnldGnDB(sGns,sGns3UTRs,sLncrnas,fOutPCGenes,fOutPCGenesSbst, \
fOutUtr3cdg,fOutUtr3cdgSbst,fOutLncrna)
#--------------------------------------------------------
#~ Output:
#~ Number of protein coding transcript in autosomes and X chromosomes 
#~ obtained was: 81579, from 20100 genes
#~ Number of protein coding transcript in autosomes and X chromosomes 
#~ obtained was: 81579, from 20100 genes
#~ Number of 3'UTRs from transcripts in autosomes and X chromosomes 
#~ obtained was: 59033, from 19351 genes
#~ Number of 3'UTRs from transcripts in autosomes and X chromosomes 
#~ obtained was: 59033, from 19351 genes
#~ Number of lncRNAs from transcripts in autosomes and X chromosomes 
#~ obtained was: 23790, from 13800 genes
#~ Number of protein coding transcript in autosomes and X chromosomes 
#~ obtained was: 81579, from 20100 genes
#~ Number of protein coding transcript in autosomes and X chromosomes 
#~ obtained was: 81579, from 20100 genes
#~ Number of 3'UTRs from transcripts in autosomes and X chromosomes 
#~ obtained was: 59033, from 19351 genes
#~ Number of 3'UTRs from transcripts in autosomes and X chromosomes 
#~ obtained was: 59033, from 19351 genes
#~ Number of lncRNAs from transcripts in autosomes and X chromosomes 
#~ obtained was: 23790, from 13800 genes
#--------------------------------------------------------


########################################################
#~ Make gene object from a fasta file in ENSEMBL format and create
# sql databases.
########################################################
#--------------------------------------------------------
#~ Set miRBase-formatted miRNA data
fstasMirnaBsFl = os.path.join(dbFldr,'mature_hg19.fa')#miRBase hg19 
# database 
srtdMirnaFms = sorted(['-'.join(l.split()[0].split('-')[1:]) for l \
in open(fstasMirnaBsFl).read().split('>') if l.strip()])
#--------------------------------------------------------
#~ All protein coding isoforms -sql file-
fOutPCTrscrptsSql = os.path.join(dbFldr,'hg19_ENSEMBL_AllPCT_Isfrms.db')
#~ Subset of protein coding isoforms -All annotated by GENCODE- -sql file-
fOutPCTrscrptsSbstSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODEPCT_Isfrms.db')
#~ All 3'UTRs from protein coding isoforms -sql file-
fOutUtr3cdgTrnscrptsSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_All3UTR_Isfrms.db')
#~ Subset of 3'UTRs from protein coding isoforms -sql file-
fOutUtr3cdgTrnscrptsSbstSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODE3UTR_Isfrms.db')
#~ Subset of lncRNA isoforms -All annotated by GENCODE- -sql file-
fOutLncrnaTrnscrptsSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODELncrns_Isfrms.db')
#~ All protein coding genes -sql file-
fOutPCGenesSql = os.path.join(dbFldr,'hg19_ENSEMBL_AllPCT_Genes.db')
#~ Subset of protein coding genes -sql file-
fOutPCGenesSbstSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODEPCT_Genes.db')
#~ All 3'UTRs from protein coding genes -sql file-
fOutUtr3cdgSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_All3UTR_Genes.db')
#~ Subset of 3'UTRs from protein coding genes -sql file-
fOutUtr3cdgSbstSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODE3UTR_Genes.db')
#~ Subset of lncRNA genes -All annotated by GENCODE- -sql file-
fOutLncrnaSql = os.path.join(dbFldr, \
'hg19_ENSEMBL_GENCODELncrns_Genes.db')
#--------------------------------------------------------
#~ Build objects for protein coding isoforms -also sql file-
lObjctsPCTrscrpts = __object__.mkArryGnObjcts(os.path.split \
(fOutPCTrscrpts)[1],srtdMirnaFms,dbFldr,sqlFl=fOutPCTrscrptsSql)
#~ Build objects of protein coding isoforms -All annotated by GENCODE- 
# -also sql file-
lObjctsPCTrscrptsSbst = __object__.mkArryGnObjcts(os.path.split \
(fOutPCTrscrptsSbst)[1],srtdMirnaFms,dbFldr,sqlFl=fOutPCTrscrptsSbstSql)
#~ Build objects for all 3'UTRs from protein coding isoforms -also sql 
# file-
lObjctsUtr3cdgTrnscrpts = __object__.mkArryGnObjcts(os.path.split \
(fOutUtr3cdgTrnscrpts)[1],srtdMirnaFms,dbFldr, \
sqlFl=fOutUtr3cdgTrnscrptsSql)
#~ Build objects for a subset of 3'UTRs from protein coding isoforms 
# -also sql file-
lObjctsUtr3cdgTrnscrptsSbst = __object__.mkArryGnObjcts( \
os.path.split(fOutUtr3cdgTrnscrptsSbst)[1],srtdMirnaFms,dbFldr, \
sqlFl=fOutUtr3cdgTrnscrptsSbstSql)
#~ Build objects for a subset of lncRNA isoforms -All annotated by 
# GENCODE- -sql file-
lObjctsLncrnaTrnscrpts = __object__.mkArryGnObjcts(os.path.split \
(fOutLncrnaTrnscrpts)[1],srtdMirnaFms,dbFldr, \
sqlFl=fOutLncrnaTrnscrptsSql)
#~ Build objects for all protein coding genes -sql file-
lObjctsPCGenes = __object__.mkArryGnObjcts(os.path.split(fOutPCGenes)[1], \
srtdMirnaFms,dbFldr,sqlFl=fOutPCGenesSql)
#~ Build objects for a subset of protein coding genes -sql file-
lObjctsPCGenesSbst = __object__.mkArryGnObjcts(os.path.split \
(fOutPCGenesSbst)[1],srtdMirnaFms,dbFldr,sqlFl=fOutPCGenesSbstSql)
#~ Build objects for all 3'UTRs from protein coding genes -sql file-
lObjctsUtr3cdg = __object__.mkArryGnObjcts(os.path.split(fOutUtr3cdg)[1], \
srtdMirnaFms,dbFldr,sqlFl=fOutUtr3cdgSql)
#~ Build objects for a subset of 3'UTRs from protein coding genes -sql 
# file-
lObjctsUtr3cdgSbst = __object__.mkArryGnObjcts(os.path.split \
(fOutUtr3cdgSbst)[1],srtdMirnaFms,dbFldr,sqlFl=fOutUtr3cdgSbstSql)
#~ Build objects for a subset of lncRNA genes -All annotated by GENCODE- 
# -sql file-
lObjctsLncrna = __object__.mkArryGnObjcts(os.path.split(fOutLncrna)[1], \
srtdMirnaFms,dbFldr,sqlFl=fOutLncrnaSql)


########################################################
#~ Make gene object from a fasta file in ENSEMBL format and test 
# congruency -with no sql database creation-.
########################################################
#~ Build objects for all 3'UTRs from protein coding isoforms 
lObjctsUtr3cdgTrnscrpts1 = __object__.mkArryGnObjcts(fOutUtr3cdgTrnscrpts, \
srtdMirnaFms,dbFldr)
lAttrbts = ['chr', 'cmmNm', 'len', 'name', 'pos', 'strnd']
for attrb in lAttrbts:
	assert [getattr(gn,attrb) for gn in lObjctsUtr3cdgTrnscrpts1] == \
	[getattr(gn,attrb) for gn in lObjctsUtr3cdgTrnscrpts]


########################################################
#~ Make gene object from another gene object (obtain subset) and test 
# congruency.
########################################################
sSubSt = set([gn.name for gn in lObjctsPCGenesSbst])
lObjctsPCGenesSbst1 = __object__.mkArrySubStObjcts(lObjctsPCGenes,sSubSt)
lAttrbts = ['chr', 'cmmNm', 'len', 'name', 'pos', 'strnd']
for attrb in lAttrbts:
	assert [getattr(gn,attrb) for gn in lObjctsPCGenesSbst1] == \
	[getattr(gn,attrb) for gn in lObjctsPCGenesSbst]

########################################################
#~ Make gene object from an sqlite3 database (obtain subset).
########################################################
sGnIntrst = set([gn.name for gn in lObjctsUtr3cdgSbst])
lObjctsUtr3cdgSbst1 = __object__.mkArrySqlObjcts(fOutUtr3cdgSql, \
srtdMirnaFms,sGnIntrst)
lAttrbts = ['chr', 'cmmNm', 'len', 'name', 'pos', 'strnd']
for attrb in lAttrbts:
	assert [getattr(gn,attrb) for gn in lObjctsUtr3cdgSbst1] == \
	[getattr(gn,attrb) for gn in lObjctsUtr3cdgSbst]

 
########################################################
#~ Make gene object from a table and test congruency.
########################################################
#--------------------------------------------------------
#~ Create a dummy tsv table
fTblTst = os.path.join(dbFldr,'hg19_ENSEMBL_GENCODELncrns_Isfrms.tsv')
tblPCGenes = []
for echL in open(fOutLncrnaTrnscrpts,'r'):
	if echL.strip() and echL[0]=='>':
		echL = echL[1:].splitlines()[0].split('|')
		gnNm,othrNm,cmmNm,chr,strts,ends,strnd = echL
		tblPCGenes.append('\t'.join([ends,strnd,gnNm,cmmNm,strts,chr]))
#
fTblTst = open(fTblTst,'w')	
fTblTst.write('\n'.join(tblPCGenes))
fTblTst.close()
#--------------------------------------------------------
#~ Make an object based on the table
fTblTst = 'hg19_ENSEMBL_GENCODELncrns_Isfrms.tsv'
lObjctsLncrnaTrnscrpts1 = __object__.mkArryGnObjcts(fTblTst, \
srtdMirnaFms,dbFldr,dType='tsv',char='\t',posNm=2,posCmmnNm=3, \
posChr=5,posStrt=4,posEnd=0,posStrnd=1)
#--------------------------------------------------------
#~ Test for congruency
lAttrbts = ['chr', 'cmmNm', 'len', 'name', 'pos', 'strnd']
for attrb in lAttrbts:
	assert [getattr(gn,attrb) for gn in lObjctsLncrnaTrnscrpts1] == \
	[getattr(gn,attrb) for gn in lObjctsLncrnaTrnscrpts]



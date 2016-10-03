#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  __db__.py part of cellery (ceRNAs linking inference)
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
Download or obtain databases.
"""

########################################################
#~ Import libraries.
########################################################
import urllib2,urllib
import gzip


########################################################
#~ Set default variables
########################################################
"""
The following variables are part of the xml code to retrieve sequences 
from ENSEMBL Biomart. biomarturlHg19 is the code for the biomart 
version. xmlElmntsHdrDlft is the code for the header. xmlElmntsBttmDlft 
is the code for the bottom. xmlElmntsDB is the code for the database."
"""
biomarturlHg19 = 'http://grch37.ensembl.org/biomart/martservice?'
xmlElmntsHdrDlft = '<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'
xmlElmntsBttmDlft = '</Dataset></Query>'


########################################################
#~ Retrive from ENSEMBL database 3'UTRs from all protein coding genes.
########################################################
def rtrn3UTRGnsHg19(fOutUtr3cdg,xmlElmntsHdr=xmlElmntsHdrDlft, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsBiotype=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False):
	"""
	Input: fOutUtr3cdg is the full path of file to save 3'UTR sequences 
	of all protein coding genes present in ENSEMBL. All of the following 
	variables are part of the xml code to retrieve sequences from 
	ENSEMBL Biomart. biomarturlHg19 is the code for the biomart version. 
	xmlElmntsHdrDlft is the code for the header. xmlElmntsBttmDlft is 
	the code for the bottom. xmlElmntsDB is the code for the database. 
	xmlElmntsBiotype is the code for the biotype. xmlElmntsLChrs is the 
	code for the list of chromosomes. xmlElmntsGns is the code for the 
	genes. xmlElmntsTrnscpts is the code for the transcripts. 
	xmlElmntsDtype is the code for the data type (i.e. 3UTR or cdna). 
	xmlElmntsExtrnlNm is the code for the header. xmlElmntsChrs is the 
	code for the chromosomes. xmlElmntsStrts  is the code for the 
	starting positions. xmlElmntsEnds is the code for the ending 
	positions. xmlElmntsStrnd is the code for the strand.
	Ouput: Fasta formated file to save 3'UTR sequences of coding genes.
	And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "3utr" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "3_utr_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "3_utr_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
	xmlElmntsLChrs,xmlElmntsGns,xmlElmntsTrnscpts,xmlElmntsDtype, \
	xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds, \
	xmlElmntsStrnd,xmlElmntsBttm])
	data = urllib.urlencode({'query': xmlElmnts})
	response = urllib2.urlopen(biomarturl, data)
	seq = response.read()
	salelseq = ['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
	splitlines()[1:])) for fs in seq.split('>') if fs.strip() and fs. \
	find('annotated')==-1 and fs.find('unavailable')==-1]
	utr3cdg = open(fOutUtr3cdg,'w')
	utr3cdg.write('\n'.join(sorted(salelseq)))
	utr3cdg.close()
	sGenes = set([x.split('>')[1].split('|')[0] for x in salelseq if \
	x.strip()])
	print "Number of 3'UTRs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database 3'UTRs from a set of protein coding 
# genes of interest.
########################################################
def rtrn3UTRGnsHg19SubSt(sGns3UTRs,fOutUtr3cdg, \
	xmlElmntsHdr=xmlElmntsHdrDlft,xmlElmntsBttm=xmlElmntsBttmDlft, \
	biomarturl=biomarturlHg19,xmlElmntsDB=False,xmlElmntsBiotype=False, \
	xmlElmntsLChrs=False,xmlElmntsGns=False,xmlElmntsTrnscpts=False, \
	xmlElmntsDtype=False,xmlElmntsExtrnlNm=False,xmlElmntsChrs=False, \
	xmlElmntsStrts=False,xmlElmntsEnds=False,xmlElmntsStrnd=False, \
	xmlElmntsLGns=False):
	"""
	Input: sGns3UTRs is the set of genes of interest. fOutUtr3cdg is the 
	full path of file to save 3'UTR sequences of protein coding genes 
	present as in sGns3UTRs. All of the following variables are part of 
	the xml code to retrieve sequences from ENSEMBL Biomart. 
	biomarturlHg19 is the code for the biomart version. xmlElmntsHdrDlft 
	is the code for the header. xmlElmntsBttmDlft is the code for the 
	bottom. xmlElmntsDB is the code for the database. xmlElmntsBiotype 
	is the code for the biotype. xmlElmntsLChrs is the code for the list 
	of chromosomes. xmlElmntsGns is the code for the genes. 
	xmlElmntsTrnscpts is the code for the transcripts. xmlElmntsDtype is 
	the code for the data type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is 
	the code for the header. xmlElmntsChrs is the code for the 
	chromosomes. xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand. xmlElmntsLGns is the code for the list 
	of genes/transcripts.
	Ouput: Fasta formated file to save 3'UTR sequences of prtoein coding 
	genes in sGns3UTRs. And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_gene_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "3utr" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "3_utr_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "3_utr_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	srtdSGnsExprs = sorted(sGns3UTRs)
	salelseq=[]
	for cnt in range(0,len(srtdSGnsExprs),400):
		lGns=srtdSGnsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsGns, \
		xmlElmntsTrnscpts,xmlElmntsDtype,xmlElmntsExtrnlNm, \
		xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds,xmlElmntsStrnd, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		seq=response.read()
		frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	utr3cdg=open(fOutUtr3cdg,'w')
	utr3cdg.write('\n'.join(sorted(salelseq)))
	utr3cdg.close()
	sGenes = set([x.split('>')[1].split('|')[0] for x in salelseq if \
	x.strip()])
	print "Number of 3'UTRs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database 3'UTRs from all protein coding 
# isoforms.
########################################################
def rtrn3UTRhg19Trnscrpts(fOutUtr3cdgTrnscrpts, \
	xmlElmntsHdr=xmlElmntsHdrDlft,xmlElmntsBttm=xmlElmntsBttmDlft, \
	biomarturl=biomarturlHg19,xmlElmntsDB=False,xmlElmntsBiotype=False, \
	xmlElmntsLChrs=False,xmlElmntsGns=False,xmlElmntsTrnscpts=False, \
	xmlElmntsDtype=False,xmlElmntsExtrnlNm=False,xmlElmntsChrs=False, \
	xmlElmntsStrts=False,xmlElmntsEnds=False,xmlElmntsStrnd=False):
	"""
	Input: fOutUtr3cdgTrnscrpts is the full path of file to save 3'UTR 
	of all protein coding isoforms present in ENSEMBL. All of the 
	following variables are part of the xml code to retrieve sequences 
	from ENSEMBL Biomart. biomarturlHg19 is the code for the biomart 
	version. xmlElmntsHdrDlft is the code for the header. 
	xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is the 
	code for the database. xmlElmntsBiotype is the code for the biotype. 
	xmlElmntsLChrs is the code for the list of chromosomes. xmlElmntsGns 
	is the code for the genes. xmlElmntsTrnscpts is the code for the 
	transcripts. xmlElmntsDtype is the code for the data type (i.e. 3UTR 
	or cdna). xmlElmntsExtrnlNm is the code for the header. 
	xmlElmntsChrs is the code for the chromosomes. xmlElmntsStrts is the 
	code for the starting positions. xmlElmntsEnds is the code for the 
	ending positions. xmlElmntsStrnd is the code for the strand.
	Ouput: Fasta formated file to save 3'UTR sequences of all protein 
	coding isoforms. And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "3utr" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "3_utr_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "3_utr_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
	xmlElmntsLChrs,xmlElmntsTrnscpts,xmlElmntsGns,xmlElmntsDtype, \
	xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds, \
	xmlElmntsStrnd,xmlElmntsBttm])
	data=urllib.urlencode({'query': xmlElmnts})
	response=urllib2.urlopen(biomarturl, data)
	seq=response.read()
	salelseq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
	splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
	fs.find('annotated')==-1 and fs.find('unavailable')==-1]
	utr3cdg=open(fOutUtr3cdgTrnscrpts,'w')
	utr3cdg.write('\n'.join(sorted(salelseq)))
	utr3cdg.close()
	sGenes = set([x.split('>')[1].split('|')[1] for x in salelseq if \
	x.strip()])
	print "Number of 3'UTRs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database 3'UTRs from a set of protein coding 
# isoforms of interest.
########################################################
def rtrn3UTRhg19TrnscrptsSubSt(sTrnscrpts3UTRs, \
	fOutUtr3cdgTrnscrpts,xmlElmntsHdr=xmlElmntsHdrDlft, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsBiotype=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False,xmlElmntsLGns=False):
	"""
	Input: sTrnscrpts3UTRs is the set of isoforms of interest. 
	fOutUtr3cdgTrnscrpts is the full path of file to save 3'UTR 
	sequences of protein coding isoforms in sTrnscrpts3UTRs. All of the 
	following variables are part of the xml code to retrieve sequences 
	from ENSEMBL Biomart. biomarturlHg19 is the code for the biomart 
	version. xmlElmntsHdrDlft is the code for the header. 
	xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is the 
	code for the database. xmlElmntsBiotype is the code for the biotype. 
	xmlElmntsLChrs is the code for the list of chromosomes. 
	xmlElmntsGns is the code for the genes. xmlElmntsTrnscpts is the 
	code for the transcripts. xmlElmntsDtype is the code for the data 
	type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is the code for the 
	header. xmlElmntsChrs is the code for the chromosomes. 
	xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand. xmlElmntsLGns is the code for the list 
	of genes/transcripts.
	Ouput: Fasta formated file to save 3'UTR sequences of coding 
	isoforms with sequences in sTrnscrpts3UTRs. And prints number of 
	transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_transcript_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "3utr" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "3_utr_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "3_utr_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	srtdsTrnscrptsExprs = sorted(sTrnscrpts3UTRs)
	salelseq=[]
	for cnt in range(0,len(srtdsTrnscrptsExprs),400):
		lGns=srtdsTrnscrptsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsTrnscpts, \
		xmlElmntsGns,xmlElmntsDtype,xmlElmntsExtrnlNm, \
		xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds,xmlElmntsStrnd, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		seq=response.read()
		frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	utr3cdg=open(fOutUtr3cdgTrnscrpts,'w')
	utr3cdg.write('\n'.join(sorted(salelseq)))
	utr3cdg.close()
	sGenes = set([x.split('>')[1].split('|')[1] for x in salelseq if \
	x.strip()])
	print "Number of 3'UTRs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database lncRNAs from a set of protein coding 
# genes of interest.
########################################################
def rtrnLncRNAGnsSubSt(sLncrnas,fOutLncrna,xmlElmntsHdr= \
	xmlElmntsHdrDlft,xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl= \
	biomarturlHg19,xmlElmntsDB=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False,xmlElmntsLGns=False):
	"""
	Input: sLncrnas is the set of lncRNA genes of interest. fOutLncrna 
	is the full path of the file to save the lncRNA gene sequences. All 
	of the following variables are part of the xml code to retrieve 
	sequences from ENSEMBL Biomart. biomarturlHg19 is the code for the 
	biomart version. xmlElmntsHdrDlft is the code for the header. 
	xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is the 
	code for the database. xmlElmntsBiotype is the code for the biotype. 
	xmlElmntsLChrs is the code for the list of chromosomes. xmlElmntsGns 
	is the code for the genes. xmlElmntsTrnscpts is the code for the 
	transcripts. xmlElmntsDtype is the code for the data type (i.e. 3UTR 
	or cdna). xmlElmntsExtrnlNm is the code for the header. 
	xmlElmntsChrs is the code for the chromosomes. xmlElmntsStrts is the 
	code for the starting positions. xmlElmntsEnds is the code for the 
	ending positions. xmlElmntsStrnd is the code for the strand. 
	xmlElmntsLGns is the code for the list of genes/transcripts.
	Ouput: Fasta formated file with lncRNA gene sequences in sLncrnas. 
	And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_gene_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	srtdLncrnas = sorted(sLncrnas)
	salelseq=[]
	for cnt in range(0,len(srtdLncrnas),400):
		lGns = srtdLncrnas[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsLGns% \
		','.join(lGns),xmlElmntsLChrs,xmlElmntsGns,xmlElmntsTrnscpts, \
		xmlElmntsDtype,xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts, \
		xmlElmntsEnds,xmlElmntsStrnd,xmlElmntsBttm])
		data = urllib.urlencode({'query': xmlElmnts})
		response = urllib2.urlopen(biomarturl, data)
		seq = response.read()
		frmtdSeq = ['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	salef = open(fOutLncrna,'w')
	salef.write('\n'.join(salelseq))
	salef.close()
	sGenes = set([x.split('>')[1].split('|')[0] for x in salelseq if \
	x.strip()])
	print "Number of lncRNAs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database lncRNAs from a set of protein coding 
# isoforms of interest.
########################################################
def rtrnLncRNATrnscrptsSubSt(sLncrnas,fOutLncrnaTrnscrpts,xmlElmntsHdr= \
	xmlElmntsHdrDlft,xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl= \
	biomarturlHg19,xmlElmntsDB=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False,xmlElmntsLGns=False):
	"""
	Input: sLncrnas is the set of lncRNA isoforms of interest. 
	fOutLncrnaTrnscrpts is the name of file to save lncRNA isoform 
	sequences. All of the following variables are part of the xml code 
	to retrieve sequences from ENSEMBL Biomart. biomarturlHg19 is the 
	code for the biomart version. xmlElmntsHdrDlft is the code for the 
	header. xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is 
	the code for the database. xmlElmntsBiotype is the code for the 
	biotype. xmlElmntsLChrs is the code for the list of chromosomes. 
	xmlElmntsGns is the code for the genes. xmlElmntsTrnscpts is the 
	code for the transcripts. xmlElmntsDtype is the code for the data 
	type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is the code for the 
	header. xmlElmntsChrs is the code for the chromosomes. 
	xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand. xmlElmntsLGns is the code for the list 
	of genes/transcripts.
	Ouput: Fasta formated file to save lncRNA sequences. And prints number 
	of transcripts and genes.
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_transcript_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	srtdLncrnasTrnscrpt = sorted(sLncrnas)
	salelseq=[]
	for cnt in range(0,len(srtdLncrnasTrnscrpt),400):
		lGns=srtdLncrnasTrnscrpt[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsLGns% \
		','.join(lGns),xmlElmntsLChrs,xmlElmntsTrnscpts,xmlElmntsGns, \
		xmlElmntsDtype,xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts, \
		xmlElmntsEnds,xmlElmntsStrnd,xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		seq=response.read()
		frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	salef=open(fOutLncrnaTrnscrpts,'w')
	salef.write('\n'.join(salelseq))
	salef.close()
	sGenes = set([x.split('>')[1].split('|')[1] for x in salelseq if \
	x.strip()])
	print "Number of lncRNAs from transcripts in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database for all protein coding genes.
########################################################
def rtrnPCGenes(fOutPCGenes,xmlElmntsHdr=xmlElmntsHdrDlft, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsBiotype=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False):
	"""
	Input: fOutPCGenes is the name of file to save protein coding gene 
	sequences. All of the following variables are part of the xml code 
	to retrieve sequences from ENSEMBL Biomart. biomarturlHg19 is the 
	code for the biomart version. xmlElmntsHdrDlft is the code for the 
	header. xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is 
	the code for the database. xmlElmntsBiotype is the code for the 
	biotype. xmlElmntsLChrs is the code for the list of chromosomes. 
	xmlElmntsGns is the code for the genes. xmlElmntsTrnscpts is the 
	code for the transcripts. xmlElmntsDtype is the code for the data 
	type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is the code for the 
	header. xmlElmntsChrs is the code for the chromosomes. 
	xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand.
	Ouput: Fasta formated file to save protein coding gene sequences. 
	And prints number of transcripts and genes.
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'		
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	salelseq=[]
	xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
	xmlElmntsLChrs,xmlElmntsGns,xmlElmntsTrnscpts,xmlElmntsDtype, \
	xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds, \
	xmlElmntsStrnd,xmlElmntsBttm])
	data=urllib.urlencode({'query': xmlElmnts})
	response=urllib2.urlopen(biomarturl, data)
	seq=response.read()
	frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
	splitlines()[1:])) for fs in seq.split('>') if fs.strip() and fs. \
	find('annotated')==-1 and fs.find('unavailable')==-1]
	salelseq.extend(frmtdSeq)
	salef=open(fOutPCGenes,'w')
	salef.write('\n'.join(salelseq))
	salef.close()
	sGenes = set([x.split('>')[1].split('|')[0] for x in salelseq if \
	x.strip()])
	print "Number of protein coding transcript in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database from a set of protein coding genes.
########################################################
def rtrnPCGenesSubSt(sGns,fOutPCGenes,xmlElmntsHdr=xmlElmntsHdrDlft, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsBiotype=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False,xmlElmntsLGns=False):
	"""
	Input: sGns is the set of genes of interest. fOutPCGenes is the name 
	of file to save protein coding gene sequences. All of the following 
	variables are part of the xml code to retrieve sequences from 
	ENSEMBL Biomart. biomarturlHg19 is the code for the biomart version. 
	xmlElmntsHdrDlft is the code for the header. xmlElmntsBttmDlft is 
	the code for the bottom. xmlElmntsDB is the code for the database. 
	xmlElmntsBiotype is the code for the biotype. xmlElmntsLChrs is the 
	code for the list of chromosomes. xmlElmntsGns is the code for the 
	genes. xmlElmntsTrnscpts is the code for the transcripts. 
	xmlElmntsDtype is the code for the data type (i.e. 3UTR or cdna). 
	xmlElmntsExtrnlNm is the code for the header. xmlElmntsChrs is the 
	code for the chromosomes. xmlElmntsStrts  is the code for the 
	starting positions. xmlElmntsEnds is the code for the ending 
	positions. xmlElmntsStrnd is the code for the strand. xmlElmntsLGns 
	is the code for the list of genes/transcripts.
	Ouput: Fasta formated file to save protein coding gene sequences. 
	And prints number of transcripts and genes.
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_gene_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	srtdSGnsExprs = sorted(sGns)
	salelseq=[]
	for cnt in range(0,len(srtdSGnsExprs),400):
		lGns=srtdSGnsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsGns, \
		xmlElmntsTrnscpts,xmlElmntsDtype,xmlElmntsExtrnlNm, \
		xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds,xmlElmntsStrnd, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		seq=response.read()
		frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	salef=open(fOutPCGenes,'w')
	salef.write('\n'.join(salelseq))
	salef.close()
	sGenes = set([x.split('>')[1].split('|')[0] for x in salelseq if \
	x.strip()])
	print "Number of protein coding transcript in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database for all protein coding isoforms.
########################################################
def rtrnPCTrscrpts(fOutPCTrscrpts,xmlElmntsHdr=xmlElmntsHdrDlft, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsBiotype=False,xmlElmntsLChrs=False, \
	xmlElmntsGns=False,xmlElmntsTrnscpts=False,xmlElmntsDtype=False, \
	xmlElmntsExtrnlNm=False,xmlElmntsChrs=False,xmlElmntsStrts=False, \
	xmlElmntsEnds=False,xmlElmntsStrnd=False):
	"""
	Input: fOutPCTrscrpts is the name of file to save protein coding 
	isoform sequences. All of the following variables are part of the 
	xml code to retrieve sequences from ENSEMBL Biomart. biomarturlHg19 
	is the code for the biomart version. xmlElmntsHdrDlft is the code 
	for the header. xmlElmntsBttmDlft is the code for the bottom. 
	xmlElmntsDB is the code for the database. xmlElmntsBiotype is the 
	code for the biotype. xmlElmntsLChrs is the code for the list of 
	chromosomes. xmlElmntsGns is the code for the genes. 
	xmlElmntsTrnscpts is the code for the transcripts. xmlElmntsDtype is 
	the code for the data type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is 
	the code for the header. xmlElmntsChrs is the code for the 
	chromosomes. xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand.
	Ouput: Fasta formated file to save protein coding transcript 
	sequences. And prints number of transcripts and genes.
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'		
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	salelseq=[]
	xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
	xmlElmntsLChrs,xmlElmntsTrnscpts,xmlElmntsGns,xmlElmntsDtype, \
	xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds, \
	xmlElmntsStrnd,xmlElmntsBttm])
	data=urllib.urlencode({'query': xmlElmnts})
	response=urllib2.urlopen(biomarturl, data)
	seq=response.read()
	frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
	splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
	fs.find('annotated')==-1 and fs.find('unavailable')==-1]
	salelseq.extend(frmtdSeq)
	salef=open(fOutPCTrscrpts,'w')
	salef.write('\n'.join(salelseq))
	salef.close()
	sGenes = set([x.split('>')[1].split('|')[1] for x in salelseq if \
	x.strip()])
	print "Number of protein coding transcript in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive from ENSEMBL database from a set of protein coding isoforms.
########################################################
def rtrnPCTrscrptsSubSt(sTrnscrpts,fOutPCTrscrpts,xmlElmntsHdr= \
	xmlElmntsHdrDlft,xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl= \
	biomarturlHg19,xmlElmntsDB=False,xmlElmntsBiotype=False, \
	xmlElmntsLChrs=False,xmlElmntsGns=False,xmlElmntsTrnscpts=False, \
	xmlElmntsDtype=False,xmlElmntsExtrnlNm=False,xmlElmntsChrs=False, \
	xmlElmntsStrts=False,xmlElmntsEnds=False,xmlElmntsStrnd=False, \
	xmlElmntsLGns=False):
	"""
	Input: sTrnscrpts is the set of isoforms of interest. fOutPCTrscrpts 
	is the full path of file to save sequences of protein coding 
	isoforms in sTrnscrpts. All of the following variables are part of 
	the xml code to retrieve sequences from ENSEMBL Biomart. 
	biomarturlHg19 is the code for the biomart version. xmlElmntsHdrDlft 
	is the code for the header. xmlElmntsBttmDlft is the code for the 
	bottom. xmlElmntsDB is the code for the database. xmlElmntsBiotype 
	is the code for the biotype. xmlElmntsLChrs is the code for the list 
	of chromosomes. xmlElmntsGns is the code for the genes. 
	xmlElmntsTrnscpts is the code for the transcripts. xmlElmntsDtype is 
	the code for the data type (i.e. 3UTR or cdna). xmlElmntsExtrnlNm is 
	the code for the header. xmlElmntsChrs is the code for the 
	chromosomes. xmlElmntsStrts  is the code for the starting positions. 
	xmlElmntsEnds is the code for the ending positions. xmlElmntsStrnd 
	is the code for the strand. xmlElmntsLGns is the code for the list 
	of genes/transcripts.
	Ouput: Fasta formated file to save sequences of coding isoforms with 
	sequences in sTrnscrpts. And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsBiotype:
		xmlElmntsBiotype = '<Filter name = "transcript_biotype" value = "protein_coding"/>'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_transcript_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsDtype:
		xmlElmntsDtype = '<Attribute name = "cdna" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsStrts:
		xmlElmntsStrts = '<Attribute name = "exon_chrom_start" />'
	if not xmlElmntsEnds:
		xmlElmntsEnds = '<Attribute name = "exon_chrom_end" />'
	if not xmlElmntsStrnd:
		xmlElmntsStrnd = '<Attribute name = "strand" />'
	#~ 
	srtdsTrnscrptsExprs = sorted(sTrnscrpts)
	salelseq=[]
	for cnt in range(0,len(srtdsTrnscrptsExprs),400):
		lGns=srtdsTrnscrptsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB,xmlElmntsBiotype, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsTrnscpts, \
		xmlElmntsGns,xmlElmntsDtype,xmlElmntsExtrnlNm, \
		xmlElmntsChrs,xmlElmntsStrts,xmlElmntsEnds,xmlElmntsStrnd, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		seq=response.read()
		frmtdSeq=['>%s\n%s'%(fs.splitlines()[0],''.join(fs. \
		splitlines()[1:])) for fs in seq.split('>') if fs.strip() and \
		fs.find('annotated')==-1 and fs.find('unavailable')==-1]
		salelseq.extend(frmtdSeq)
	cdg=open(fOutPCTrscrpts,'w')
	cdg.write('\n'.join(sorted(salelseq)))
	cdg.close()
	sGenes = set([x.split('>')[1].split('|')[1] for x in salelseq if \
	x.strip()])
	print "Number of protein coding transcript in autosomes and X chromosomes obtained was: %s, from %s genes"% \
	(len(salelseq),len(sGenes))
	return 0


########################################################
#~ Retrive annotations from ENSEMBL database for a set of genes.
########################################################
def rtrnAnntGenesSubSt(sGns,fOutAnntGns,xmlElmntsHdr=False, \
	xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsLChrs=False,xmlElmntsGns=False, \
	xmlElmntsTrnscpts=False,xmlElmntsExtrnlNm=False,xmlElmntsChrs=False, \
	xmlElmntsLGns=False,xmlElmntsGO=False):
	"""
	Input: sGns is the set of isoforms of interest. fOutAnntGns is the 
	full path of '.tsv' file to save annotations of the genes in sGns. 
	All of the following variables are part of the xml code to retrieve 
	sequences from ENSEMBL Biomart. biomarturlHg19 is the code for the 
	biomart version. xmlElmntsHdrDlft is the code for the header. 
	xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is the 
	code for the database. xmlElmntsLChrs is the code for the list of 
	chromosomes. xmlElmntsGns is the code for the genes. 
	xmlElmntsTrnscpts is the code for the transcripts. xmlElmntsExtrnlNm 
	is the code for the header. xmlElmntsChrs is the code for the 
	chromosomes. xmlElmntsLGns is the code for the list 
	of genes/transcripts. xmlElmntsGO is the code for the GO terms.
	Ouput: '.tsv' formated file to save the annotations for the genes in 
	sGns. And prints number of transcripts and genes.
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsHdr:
		xmlElmntsHdr = '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_gene_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsGO:
		xmlElmntsGO = '<Attribute name = "name_1006" />'
	#--------------------------
	#~ Retrieve annotation from BioMart
	srtdsGnsExprs = sorted(sGns)
	salel=[]
	for cnt in range(0,len(srtdsGnsExprs),400):
		lGns=srtdsGnsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsGns, \
		xmlElmntsTrnscpts,xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsGO, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		annttn=response.read()
		frmtdAnnttn = [x.split('\t') for x in annttn.splitlines() if \
		x.strip()]
		salel.extend(frmtdAnnttn)
	#--------------------------
	#~ Group annotations
	dTrnscrptsAnnt = {}
	for gn,trnscrpt,extrnlNm,chr,GOtrm in salel:
		if dTrnscrptsAnnt.has_key((gn,trnscrpt,extrnlNm,chr)):
			dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr].append(GOtrm)
		else:
			dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr]=[GOtrm]
	#--------------------------
	#~ Join annotations
	salAnttn = []
	for gn,trnscrpt,extrnlNm,chr in sorted(dTrnscrptsAnnt.keys()):
		salAnttn.append('\t'.join([gn,trnscrpt,extrnlNm,chr]+['.'.join( \
		dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr])]))
	#--------------------------
	#~ Save file
	anntFlOpnd=open(fOutAnntGns,'w')
	anntFlOpnd.write('\n'.join(sorted(salAnttn)))
	anntFlOpnd.close()
	sGenes = set([x.split('\t')[0] for x in salAnttn if x.strip()])
	sGns = set([x.split('\t')[1] for x in salAnttn if x.strip()])
	print "The number of transcripts annotated in autosomes and X chromosomes was: %s from %s genes"% \
	(len(sGns),len(sGenes))
	return 0


########################################################
#~ Retrive annotations from ENSEMBL database for a set of transcripts.
########################################################
def rtrnAnntTrscrptsSubSt(sTrnscrpts,fOutAnntTrnscrpts,xmlElmntsHdr= \
	False,xmlElmntsBttm=xmlElmntsBttmDlft,biomarturl=biomarturlHg19, \
	xmlElmntsDB=False,xmlElmntsLChrs=False,xmlElmntsGns=False, \
	xmlElmntsTrnscpts=False,xmlElmntsExtrnlNm=False,xmlElmntsChrs=False, \
	xmlElmntsLGns=False,xmlElmntsGO=False):
	"""
	Input: sTrnscrpts is the set of isoforms of interest. fOutAnntTrnscrpts 
	is the full path of '.tsv' file to save annotations of isoforms in 
	sTrnscrpts. All of the following variables are part of the xml code 
	to retrieve sequences from ENSEMBL Biomart. biomarturlHg19 is the 
	code for the biomart version. xmlElmntsHdrDlft is the code for the 
	header. xmlElmntsBttmDlft is the code for the bottom. xmlElmntsDB is 
	the code for the database. xmlElmntsLChrs is the code for the list 
	of chromosomes. xmlElmntsGns is the code for the genes. 
	xmlElmntsTrnscpts is the code for the transcripts. xmlElmntsExtrnlNm 
	is the code for the header. xmlElmntsChrs is the code for the 
	chromosomes. xmlElmntsLGns is the code for the list 
	of genes/transcripts. xmlElmntsGO is the code for the GO terms.
	Ouput: '.tsv' formated file to save the annotations for the isoforms 
	in sTrnscrpts. And prints number of transcripts and genes. 
	NOTE: Retrieves sequences for grch37 == hg19 by default.
	"""
	#~ 
	if not xmlElmntsHdr:
		xmlElmntsHdr = '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'
	if not xmlElmntsDB:
		xmlElmntsDB = '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
	if not xmlElmntsLChrs:
		xmlElmntsLChrs = '<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"/>'
	if not xmlElmntsLGns:
		xmlElmntsLGns = '<Filter name = "ensembl_transcript_id" value = "%s"/>'
	if not xmlElmntsGns:
		xmlElmntsGns = '<Attribute name = "ensembl_gene_id" />'
	if not xmlElmntsTrnscpts:
		xmlElmntsTrnscpts = '<Attribute name = "ensembl_transcript_id" />'
	if not xmlElmntsExtrnlNm:
		xmlElmntsExtrnlNm = '<Attribute name = "external_gene_name" />'
	if not xmlElmntsChrs:
		xmlElmntsChrs = '<Attribute name = "chromosome_name" />'
	if not xmlElmntsGO:
		xmlElmntsGO = '<Attribute name = "name_1006" />'
	#--------------------------
	#~ Retrieve annotation from BioMart
	srtdsTrnscrptsExprs = sorted(sTrnscrpts)
	salel=[]
	for cnt in range(0,len(srtdsTrnscrptsExprs),400):
		lGns=srtdsTrnscrptsExprs[cnt:cnt+400]
		xmlElmnts = ''.join([xmlElmntsHdr,xmlElmntsDB, \
		xmlElmntsLChrs,xmlElmntsLGns%','.join(lGns),xmlElmntsTrnscpts, \
		xmlElmntsGns,xmlElmntsExtrnlNm,xmlElmntsChrs,xmlElmntsGO, \
		xmlElmntsBttm])
		data=urllib.urlencode({'query': xmlElmnts})
		response=urllib2.urlopen(biomarturl, data)
		annttn=response.read()
		frmtdAnnttn = [x.split('\t') for x in annttn.splitlines() if \
		x.strip()]
		salel.extend(frmtdAnnttn)
	#--------------------------
	#~ Group annotations
	dTrnscrptsAnnt = {}
	for gn,trnscrpt,extrnlNm,chr,GOtrm in salel:
		if dTrnscrptsAnnt.has_key((gn,trnscrpt,extrnlNm,chr)):
			dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr].append(GOtrm)
		else:
			dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr]=[GOtrm]
	#--------------------------
	#~ Join annotations
	salAnttn = []
	for gn,trnscrpt,extrnlNm,chr in sorted(dTrnscrptsAnnt.keys()):
		salAnttn.append('\t'.join([gn,trnscrpt,extrnlNm,chr]+['.'.join( \
		dTrnscrptsAnnt[gn,trnscrpt,extrnlNm,chr])]))
	#--------------------------
	#~ Save file
	anntFlOpnd=open(fAnntTrnscrpts,'w')
	anntFlOpnd.write('\n'.join(sorted(salAnttn)))
	anntFlOpnd.close()
	sGenes = set([x.split('\t')[0] for x in salAnttn if x.strip()])
	sTrnscrpts = set([x.split('\t')[1] for x in salAnttn if x.strip()])
	print "The number of transcripts annotated in autosomes and X chromosomes was: %s from %s genes"% \
	(len(sTrnscrpts),len(sGenes))
	return 0

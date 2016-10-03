#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  __object__.py part of cellery (ceRNAs linking inference)
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
Libraries to create and manipulate gene objects
"""

########################################################
#~ Import libraries.
########################################################
from cellery import __classes__
from cellery import exceptions
from numpy import empty,array,float32,float64

import copy
import os
import sqlite3
import string


########################################################
#~ Return length from an array of intervals.
########################################################
def clcLen(aMrdgIntrvls):
	"""
	Input: aMrdgIntrvls is a list of ranges in a tuple (start,end).
	Ouput: leN is the summed length of the intervals.
	"""
	leN = 0
	for strt,end in aMrdgIntrvls:
		leN += end-strt
	return leN


########################################################
#~ Merge positions of overlapping intervals for individual genes.
########################################################
def merge_ranges(ranges):
	"""
	Input: ranges is a list of ranges in a tuple (start,end).
	Ouput: current_start, current_stop is the start and end of 
	merged overlapping ranges.
	"""
	ranges = iter(sorted(ranges))
	current_start, current_stop = next(ranges)
	for start, stop in ranges:
		if start > current_stop:
			yield current_start, current_stop
			current_start, current_stop = start, stop
		else:
			current_stop = max(current_stop, stop)
	yield current_start, current_stop


########################################################
#~ Make gene object from a fasta file in ENSEMBL format or a text table.
########################################################
def mkArryGnObjcts(dtFl,srtdMirnaFms,pthFldrdtFl,dType='fas',char=False,
	posNm=0,posCmmnNm=2,posChr=3,posStrt=4,posEnd=5,posStrnd=6,
	sqlFl=False):
	"""
	Input: dtFl is the name of the input fasta file. srtdMirnaFms is 
	the sorted list of miRNA families to study. pthFldrdtFl is the 
	path to the folder with dtFl. Optionally, dType is the type of file, 
	valid types for fasta are {"fa","fasta","fas"}. char is the 
	character that splits the columns, posNm is the position of gene/
	transcript identifier/name, posCmmnNm is the position of gene/
	transcript common  name, posChr is the position of the chromosome, 
	posStrt is the position of the start values, posEnd is the position 
	of the end values, posStrnd is the position of the strands. 
	Optionally if a value for sqlFl is included a sqlite3 database will 
	be created for the input file.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas of in dtFl. Optionally if a value for sqlFl is 
	included a sqlite3 database will be created for the input file.
	"""	
	#fix parameters for process input file
	if string.lower(dType) in {"fa","fasta","fas"}:
		strtChr='>'
		if not char:
			char='|'
	else:
		strtChr = False
		try:
			assert char
		except:
			raise exceptions.CelleryWarningObjct \
			('%s file is of type TXT with column splits'%dtFl,'False')
	#process file
	dtFl = os.path.join(pthFldrdtFl,dtFl)
	#make a dictionary of merged overlapping intervals for genes
	dNmTrpleMrgdChrLIntrvlsStrndCmmnm = rtrnDNmTrpleMrgdChrLIntrvlsStrnd \
	(dtFl,char,strtChr,posNm,posCmmnNm,posChr,posStrt,posEnd,posStrnd, \
	sqlFl)
	srtdDtNms = sorted(dNmTrpleMrgdChrLIntrvlsStrndCmmnm.keys())
	cntFls = len(srtdDtNms)
	#create objects
	lDtObjcts = empty(cntFls,dtype=object)
	srtdDtNms.reverse()
	pos = -1
	cntFls -= 1
	while pos<cntFls:
		pos += 1
		name = srtdDtNms.pop()
		chr,aMrdgIntrvls,strnd,cmmNm = \
		dNmTrpleMrgdChrLIntrvlsStrndCmmnm.pop(name)
		gn = __classes__.gene(name,srtdMirnaFms)
		gn.pos = pos
		gn.chr = chr
		gn.aIntrvls = aMrdgIntrvls
		gn.len = clcLen(aMrdgIntrvls)
		gn.strnd = strnd
		gn.cmmNm = cmmNm
		lDtObjcts[pos] = gn	
	return lDtObjcts
	

########################################################
#~ Make gene object from an sqlite3 database (obtain subset).
########################################################
def mkArrySqlObjcts(sqlFl,srtdMirnaFms,sGnIntrst=False):
	"""
	Input: sqlFl is a sqlite3 database with the attributes of interest.
	srtdMirnaFms is the sorted list of miRNA families to study. 
	Optionally, sGnIntrst is a set with names of gene identifiers/names 
	of interest.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas of in sGnIntrst. If sGnIntrst is False all genes in 
	the database are going to be retrieved.
	NOTE: starts are 0-based and ends 1-based. 
	NOTE: id is the field for gene identifiers/names.
	"""	
	#make a dictionary of merged overlapping intervals for genes
	dNmTrpleMrgdChrLIntrvlsStrndCmmnm = rtrnSqlFl(sqlFl,sGnIntrst)
	srtdDtNms = sorted(dNmTrpleMrgdChrLIntrvlsStrndCmmnm.keys())
	cntFls = len(srtdDtNms)
	#create objects
	lDtObjcts = empty(cntFls,dtype=object)
	srtdDtNms.reverse()
	pos = -1
	cntFls -= 1
	while pos<cntFls:
		pos += 1
		name = srtdDtNms.pop()
		chr,aMrdgIntrvls,strnd,cmmNm = \
		dNmTrpleMrgdChrLIntrvlsStrndCmmnm.pop(name)
		gn = __classes__.gene(name,srtdMirnaFms)
		gn.pos = pos
		gn.chr = chr
		gn.aIntrvls = aMrdgIntrvls
		gn.len = clcLen(aMrdgIntrvls)
		gn.strnd = strnd
		gn.cmmNm = cmmNm
		lDtObjcts[pos] = gn	
	return lDtObjcts


########################################################
#~ Make gene object from another gene object (obtain subset).
########################################################
def mkArrySubStObjcts(lDtObjctsIn,sSubSt):
	"""
	Input: lDtObjctsIn is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas of in dtFl. sSubSt is a set of names of interest 
	contained in lDtObjctsIn and from which a new object is going to be
	build.
	Output: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas of in sSubSt.
	NOTE: All atributes for objects are copied from those in 
	lDtObjctsIn, excepting pos which is updated.
	"""	
	dGnsInlDtObjctsPos = dict([(dt.name,dt.pos) for dt in lDtObjctsIn])
	sGnsInlDtObjcts = set(dGnsInlDtObjctsPos.keys())
	sGnsDiff = sSubSt.difference(sGnsInlDtObjcts)
	sGnsIntrsct = sSubSt.intersection(sGnsInlDtObjcts)
	try:
		if sGnsDiff:
			raise exceptions.CelleryWarningObjct \
			('%s identifiers not present in'%len(sGnsDiff), \
			'for the input list of objects')
	except exceptions.CelleryWarningObjct as err:
		 print err
		 pass
	#process input set
	srtdDtNms = sorted(sGnsIntrsct)
	cntFls = len(srtdDtNms)
	#substract objects
	lDtObjcts = empty(cntFls,dtype=object)
	srtdDtNms.reverse()
	pos = -1
	cntFls -= 1
	while pos<cntFls:
		pos += 1
		name = srtdDtNms.pop()
		#deep copy gene object
		gn = copy.deepcopy(lDtObjctsIn[dGnsInlDtObjctsPos[name]])
		#update position
		gn.pos = pos
		lDtObjcts[pos] = gn	
	return lDtObjcts
	

########################################################
#~ Make a sqlite3 database for genes/lncRNAs or regions of interest.
########################################################
def mkSqlFl(dNmTrpleMrgdChrLIntrvlsStrndCmmnm,sqlFl):
	"""
	Input: dNmTrpleMrgdChrLIntrvlsStrndCmmnm is a dictionary of 
	gene/transcript identifier/name in column 0 -or posNm- as keys and 
	as values tuples of (chromosome, list of merged intervals, strand, 
	and common names). A sqlite3 database will be created for the input 
	dictionary in the file sqlFl.
	Output: A sqlite3 database will be created for the input dictionary 
	in the file sqlFl.
	NOTE: starts are 0-based and ends 1-based.
	"""
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	c.execute('''CREATE TABLE records (id TEXT, chr TEXT, strts TEXT, ends TEXT, strnd TEXT, cmmnNm TEXT)''')
	records = []
	sGns = set(dNmTrpleMrgdChrLIntrvlsStrndCmmnm.keys())
	while sGns:
		gn = sGns.pop()
		chr,aMrdgIntrvls,strnd,cmmNm = \
		dNmTrpleMrgdChrLIntrvlsStrndCmmnm.pop(gn)
		strts,ends = [],[]
		for strt,end in aMrdgIntrvls:
			strts.append(str(strt))
			ends.append(str(end))
		strts = ';'.join(strts)
		ends = ';'.join(ends)
		records.append((gn,chr,strts,ends,strnd,cmmNm))
	records = tuple(records)
	c.executemany('insert into records VALUES (?,?,?,?,?,?)', records)
	# create indexes. Decrease complexity of querying
	c.execute("CREATE INDEX index_records on records (id);")
	conn.commit()
	conn.close()
	return 0	


########################################################
#~ Make a sqlite3 database for genes/lncRNAs objects.
########################################################
def mkSqlFlMirnVls(lDtObjcts,attrbt,sqlFl):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas. attrbt is a pointer to the array with the miRNA 
	scores or counts of interest. A sqlite3 database will be created for 
	the input attribute with fields having the sames in the pointer 
	aMirnaNms.
	Output: A sqlite3 database will be created for the input dictionary 
	in the file sqlFl.
	"""
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	srtdMirnaNms = getattr(lDtObjcts[0],'aMirnaNms')
	mirnFlds = ['"%s" REAL'%mirn for mirn in srtdMirnaNms]
	c.execute('''CREATE TABLE %s (id TEXT, %s)'''% \
	(attrbt,', '.join(mirnFlds)))
	records = []
	for dt in lDtObjcts:
		mirnaAttrbt = [float64(attrbtVal) for attrbtVal in \
		getattr(dt,attrbt)]
		mirnaAttrbt.insert(0,dt.name)
		records.append(mirnaAttrbt)
	records = tuple(records)
	emptVls = ','.join(['?' for mirn in xrange(len(srtdMirnaNms)+1)])
	c.executemany('insert into %s VALUES (%s);'%(attrbt,emptVls), \
	records)
	c.execute('CREATE INDEX index_%s on %s (id);'%(attrbt,attrbt))
	conn.commit()
	conn.close()	
	return 0


########################################################
#~ Merge positions of overlapping intervals for all genes of interest.
########################################################
def mrgeDtInf(dNmTrpleSChrLIntrvlsSStrndCmmnm):
	"""
	Input: dNmTrpleSChrLIntrvlsSStrndCmmnm is a dictionary of gene names as 
	keys and as values tuples of (chromosome, list of sorted intervals, 
	and strand).
	Output: dNmTrpleMrgdChrLIntrvlsStrndCmmnm is a dictionary of gene 
	names as keys and as values tuples of (chromosome, list of merged 
	intervals, and strand).
	"""
	dNmTrpleMrgdChrLIntrvlsStrndCmmnm = {}
	sNames = set(dNmTrpleSChrLIntrvlsSStrndCmmnm.keys())
	while sNames:
		name = sNames.pop()
		sChrs,intrvls,sStrnd,sCmmNms = \
		dNmTrpleSChrLIntrvlsSStrndCmmnm.pop(name)
		try:
			if len(sChrs)!=1:
				raise exceptions.CelleryWarningObjct \
				('%s chromosomes for gene'%len(sChrs),name)
			if len(sStrnd)!=1:
				raise exceptions.CelleryWarningObjct \
				('%s strands for gene'%len(sStrnd),name)
			if len(sCmmNms)!=1:
				raise exceptions.CelleryWarningObjct \
				('%s common names for gene'%len(sCmmNms),name)
			aMrdgIntrvls = array(list(merge_ranges(intrvls)))
			chr = sChrs.pop()
			strnd = sStrnd.pop()
			cmmNm = sCmmNms.pop()
		except exceptions.CelleryWarningObjct as err:
			 print err
			 pass
		dNmTrpleMrgdChrLIntrvlsStrndCmmnm[name] = (chr,aMrdgIntrvls, \
		strnd,cmmNm)
	return dNmTrpleMrgdChrLIntrvlsStrndCmmnm


########################################################
#~ Merge overlapping intervals for genes from a fasta file in ENSEMBL 
# format or a text table.
########################################################
def rtrnDNmTrpleMrgdChrLIntrvlsStrnd(dtFl,char,strtChr,posNm,posCmmnNm, \
	posChr,posStrt,posEnd,posStrnd,sqlFl):
	"""
	Input: dtFl is an ENSEMBL fasta format or text table with the 
	following fields in the header in order = "gene/transcript 
	identifier/name, gene/transcript identifier/name, common gene/
	transcript name|starts -joined by ;-, ends -joined by ;-, strand (1 
	or -1)". char is the character that splits the columns, if strtChr 
	if included (i.e. strtChr='>') only those lines are going to be 
	scanned. posNm is the position of gene/transcript identifier/name, 
	posCmmnNm is the position of gene/transcript common  name, posChr is 
	the position of the chromosome, posStrt is the position of the start 
	values, posEnd is the position of the end values, posStrnd is the 
	position of the strands. Optionally if a value for sqlFl is included 
	a sqlite3 database will be created for the input file.
	Output: dNmTrpleMrgdChrLIntrvlsStrndCmmnm is a dictionary of 
	gene/transcript identifier/name in column 0 -or posNm- as keys and 
	as values tuples of (chromosome, list of merged intervals, and 
	strand, and common names). Optionally if a value for sqlFl is 
	included a sqlite3 database will be created for the input file.
	NOTE: input starts and ends are 1-based (as default in ENSEMBL). 
	output starts are 0-based and ends 1-based.
	"""
	dNmTrpleSChrLIntrvlsSStrndCmmnm = {}
	for eLine in open(dtFl,'r'):
		if eLine.strip():
			fullInf = eLine.strip().split(char)
			if strtChr:
				if eLine[0][0]!=strtChr:
					continue
				else:
					name = fullInf[posNm][1:]
			else:
				name = fullInf[posNm]
			cmmnNm = fullInf[posCmmnNm]
			chr = fullInf[posChr]
			starts = array(sorted([int(x) for x in \
			fullInf[posStrt].split(';')]))
			starts -= 1
			ends = array(sorted([int(x) for x in \
			fullInf[posEnd].split(';')]))
			intrvls = zip(starts,ends)
			strnd =  fullInf[posStrnd]
			if dNmTrpleSChrLIntrvlsSStrndCmmnm.has_key(name):
				dNmTrpleSChrLIntrvlsSStrndCmmnm[name][0].add(chr)
				dNmTrpleSChrLIntrvlsSStrndCmmnm[name][1].extend(intrvls)
				dNmTrpleSChrLIntrvlsSStrndCmmnm[name][2].add(strnd)
				dNmTrpleSChrLIntrvlsSStrndCmmnm[name][3].add(cmmnNm)
			else:
				dNmTrpleSChrLIntrvlsSStrndCmmnm[name] = ({chr},intrvls, \
				{strnd},{cmmnNm})	
	dNmTrpleMrgdChrLIntrvlsStrndCmmnm = \
	mrgeDtInf(dNmTrpleSChrLIntrvlsSStrndCmmnm)
	if sqlFl:
		try:
			mkSqlFl(dNmTrpleMrgdChrLIntrvlsStrndCmmnm,sqlFl)
		except:
			pass
	return dNmTrpleMrgdChrLIntrvlsStrndCmmnm


########################################################
#~ Read a sqlite3 database for genes/lncRNAs or regions of interest.
########################################################
def rtrnSqlFl(sqlFl,sGnIntrst):
	"""
	Input: sqlFl is a sqlite3 database with the attributes of interest.
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
	conn.close()
	return dNmTrpleMrgdChrLIntrvlsStrndCmmnm


########################################################
#~ Read a sqlite3 database for genes/lncRNAs objects.
########################################################
def rtrnSqlFlMirnVls(lDtObjcts,attrbt,sqlFl):
	"""
	Input: lDtObjcts is a list with gene/lncrna objects with pointers  
	"pos", "name", "aIntrvls", "chr", "cmmNm", and "strnd" for the 
	genes/lncrnas. attrbt is a pointer to the array with the miRNA 
	scores or counts of interest. sqlFl is a sqlite3 database with 
	fields having (a subset) of names in the same order of aMirnaNms.
	Output: lDtObjcts with the updated attrbt with the values stored in
	sqlFl in the same order of aMirnaNms. 
	NOTE: id is the field for gene identifiers/names.
	NOTE: Only genes information for lDtObjcts is going to be retrieved.
	"""
	srtdMirnaNms = getattr(lDtObjcts[0],'aMirnaNms')
	conn = sqlite3.connect(sqlFl)
	c = conn.cursor()
	# retrive list of miRNAs present
	sMirnsInSql = set([str(mirn[1]) for mirn in \
	c.execute("PRAGMA table_info('%s')"%attrbt) if str(mirn)!='id'])
	try:
		if not set(srtdMirnaNms).issubset(sMirnsInSql):
			sGnsDiff = set(srtdMirnaNms).difference(sMirnsInSql)
			raise exceptions.CelleryWarningObjct \
			('%s miRNAs not present in'%len(sGnsDiff),sqlFl)
	except exceptions.CelleryWarningObjct as err:
		 print err
		 pass
	# retrieve gene info
	dDtNmPos = dict([(dt.name,dt.pos) for dt in lDtObjcts])
	sDts = set(dDtNmPos.keys())
	sDtsInSql = set()
	for info in c.execute('SELECT id,%s FROM %s'% \
		(','.join(['"%s"'%m for m in srtdMirnaNms]),attrbt)):
		dtNm = str(info[0])#identifier
		if dtNm in sDts:
			dtPos = dDtNmPos[dtNm]
			aAttrbt = array([float32(val) for val in info[1:]])
			setattr(lDtObjcts[dtPos],attrbt,aAttrbt)
			sDts.remove(dtNm)		
			sDtsInSql.add(dtNm)
	# test for same genes present in input and output
	sDts = set(dDtNmPos.keys())
	try:
		sGnsDiff = sDts.difference(sDtsInSql)
		if sGnsDiff:
			raise exceptions.CelleryWarningObjct \
			('%s identifiers not present in'%len(sGnsDiff),sqlFl)
	except exceptions.CelleryWarningObjct as err:
		 print err
		 pass
	return 0

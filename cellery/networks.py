#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  networks.py part of cellery (ceRNAs linking inference)
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
Libraries to analyse and manipulate graphs and networks.
"""

########################################################
#~ Import libraries.
########################################################
from cellery import exceptions
from copy import copy
from numpy import array,float32,isnan,empty,ma,nan,percentile,vectorize

import networkx as nx
import numpy as np

#----------------------------
#~ GML based methods
#----------------------------

########################################################
#~ Annotate nodes with information of interest in GML graph files
########################################################
def anntGMLGrph(dLbldNmAnntn,G_targetGML,G_targetGMLOut=False,vrbse=True):
	"""
	Input: dLbldNmAnntn is a dictionary whose keys are annotation 
	labels/categories and their values dictionaries with isoforms/gene 
	names as keys and their annotation(s) as values. G_targetGML is the 
	graph of interest to annotate in '.gml' format. Optionally, 
	G_targetGMLOut is the output file to write the annotated graph in 
	'.gml' format. If vrbse is True, the log for the run is going to be 
	printed.
	Output: G_annttdGML (or G_targetGMLOut if True) is the output 
	annotated graph in '.gml' format with the nodes annotated following 
	dLbldNmAnntn. ovrAllogRun is the log for the run of the method.
	NOTE: If common names (i.e. values for the "cmmNm" pointer) are 
	going to be included in the graph, they shall be under the label 
	"cmmNm".
	"""
	#----------------------------
	# Method to parse the input annotations
	def prsAnnt(dNmAnnt):
		"""
		Input: dNmAnnt is a dicitonary with the name as key and the 
		annotations as values.
		Output: sAnnt is the set of all annotation included in the 
		values of dNmAnnt.
		"""
		sAnnt=set()
		for annt in dNmAnnt.values():
			sAnnt.update(annt.split('.'))
		sAnnt.remove('N')
		return sAnnt
	#----------------------------
	# Extract the set of annotations for each label
	ovrAllogRun = []
	sLblsAnntn = set(dLbldNmAnntn.keys())
	srtdLblsAnntn = sorted(sLblsAnntn)
	dLblsAnntns = {}
	sNms = set()
	for lbl in sLblsAnntn:
		dNmAnntn = dLblsAnntns[lbl]
		sNms.update(dNmAnntn.keys())
		sAnnt = prsAnnt(dNmAnnt)
		dLblsAnntns[lbl] = sorted(sAnnt)
	anntCmmnNm = False
	if 'cmmNm' in sLblsAnntn:
		anntCmmnNm = True
		tmpVrbl = dLblsAnntns.pop('cmmNm')
		del(tmpVrbl)#free memory
		try:
			raise exceptions.CelleryWarningObjct \
			('Common names are going to be annotated in', \
			'the output graph')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
	try:
		raise exceptions.CelleryWarningObjct \
		('%s annotation categories '%len(sLblsAnntn), \
		'and %s nodes are going to be included'%len(sNms))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
		pass
	#----------------------------
	# Extract isoform/gene names in the graph
	sNmsInGrph=set()
	for l in open(G_targetGML,'r'):
		if l.find('label')>-1:
			ensg=l.split('"')[1]
			sNmsInGrph.add(ensg)
	try:
		raise exceptions.CelleryWarningObjct \
		('%s'%len(sNmsInGrph),'nodes are present in the input graph')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
		pass
	assert sNmsInGrph.issubset(sNms)#names in graph are a subset
	#----------------------------
	# Include the annotations in the current graph
	salel=[]
	for l in open(G_targetGML,'r'):
		if l.find('label')>-1:
			salel.append(l)
			name = l.split('"')[1]
			#----------------------------
			# Include common name
			if anntCmmnNm:
				cmmNmVal = dLbldNmAnntn['cmmNm'].get(name,'N')
				salel.append('    name "%s"\n'%cmmNmVal)
			#----------------------------
			# Include annotations
			for lbl in srtdLblsAnntn:
				srtdAnnt = dLblsAnntns[lbl]
				sAnntInNm = dLbldNmAnntn[lbl].get(name,False)
				if sAnntInNm:
					sAnntInNm = set(sAnntInNm.split('.'))
				for annt in srtdAnnt:
					if not sAnntInNm:
						clr = 'white'
					elif annt in sAnntInNm:
						clr = 'red'
					else:
						clr = 'black'
					salel.append('    %s_%s "%s"\n'%(lbl,annt,clr))
		else:
			salel.append(l)
	#----------------------------
	# Write annotated graph
	if G_targetGMLOut:
		salef=open(G_targetGMLOut,'w')
	else:
		salef=open(G_targetGML,'w')
	salef.write(''.join(salel))
	salef.close()
	return ovrAllogRun


########################################################
#~ Annotate a subgraph made up by a subset of nodes of interest using 
# GML graph files
########################################################
def anntGMLSubGrph(sNds,G_targetGML,G_targetGMLOut,vrbse=True):
	"""
	Input: sNds is a set of nodes of interest to substract from a graph. 
	G_targetGML is the input graph of interest (from which a subset is 
	going to be obtained) in '.gml' format. G_targetGMLOut is the output 
	file to write the annotated subgraph within the supergraph in '.gml' 
	format . If vrbse is True,  the log forthe run is going to be 
	printed.
	Output: G_targetGMLOut is the output file with the annotation 
	"frstNtwrk" as 1.0 if nodes are present in sNds and 0.0 otherwise. 
	ovrAllogRun is the log of the run.
	NOTE: Assumes undirected graphs for input and output.
	"""
	#----------------------------
	#~ get nodes of interest
	sIdsIntrst = set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			if nodeNm in sNds:
				sIdsIntrst.add(idT)
	#----------------------------
	#~ set parameters
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest'%len(sIdsCnntng), \
		'are going to be retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ get gml information for nodes of interest 
	G_targetGMLOut = open(G_targetGMLOut,'w')
	G_targetGMLOut.write('graph [\n')
	hldr = []
	cntFnlNds = 0
	for l in open(G_targetGML,'r'):
		if l.find('node')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
			if idT in sIdsIntrst:
				cntFnlNds+=1
				frstNtwrk = '1.0'
			else:
				frstNtwrk = '0.0'
			hldr.append('  node [\n    id %s\n    frstNtwrk %s\n'% \
			(idT,frstNtwrk))
		elif hldr:
			hldr.append(l)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest 
	hldr = []
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst and trgt in sIdsIntrst:
				hldr.append('  edge [\n    source %s\n    target %s\n'% \
				(src,trgt))
		elif hldr:
			hldr.append(l)
	if hldr:#last case
		G_targetGMLOut.write(''.join(hldr))
	G_targetGMLOut.write(']')
	G_targetGMLOut.close()
	#----------------------------
	#~ report final values
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest'%cntFnlNds,'were annotated')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Annotate a subgraph made up by a subset of edges of interest using 
# GML graph files
########################################################
def anntGMLSubGrphEdgs(sTplSrtdEdgsFrstNtwrk,G_targetGML,G_targetGMLOut, \
	vrbse=True):
	"""
	Input: sTplSrtdEdgsFrstNtwrk is a set of sorted pair-nodes 
	representing edges of interest to substract from a graph. 
	G_targetGML is the input graph of interest (from which a subset is 
	going to be obtained) in '.gml' format. G_targetGMLOut is the output 
	file to write the annotated subset graph in '.gml' format. If vrbse 
	is True, the log for the run is going to be printed.
	Output: G_targetGMLOut is the output file with the annotation 
	"frstNtwrk" as 1.0 if edges are present in sTplSrtdEdgsFrstNtwrk and 
	0.0 otherwise. ovrAllogRun is the log of the run.
	NOTE: Assumes undirected graphs for input and output.
	"""
	#----------------------------
	#~ get nodes of in the graph
	sNds = set()
	for ndA,ndB in sTplSrtdEdgsFrstNtwrk:
		sNds.update({ndA,ndB})
	sIdsIntrst = set()
	dNodeNmIdT = {}
	sNodeNms = set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			if nodeNm in sNds:
				sIdsIntrst.add(idT)
				dNodeNmIdT[nodeNm] = idT
				sNodeNms.add(nodeNm)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get nodes of interest
	sShrdVldEdgs = []
	while len(sTplSrtdEdgsFrstNtwrk)>0:
		ndA,ndB = sTplSrtdEdgsFrstNtwrk.pop()
		if ndA in sNodeNms and ndB in sNodeNms:
			idTA,idTB = dNodeNmIdT[ndA],dNodeNmIdT[ndB]
			sShrdVldEdgs.append(sorted([idTA,idTB]))
	#----------------------------
	# report preliminary
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges and %s nodes of interest are'% \
		(len(sTplSrtdEdgsFrstNtwrk),len(sNodeNms)), \
		'going to be searched')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ get gml information for nodes of interest
	G_targetGMLOut = open(G_targetGMLOut,'w')
	G_targetGMLOut.write('graph [\n')
	hldr = []
	for l in open(G_targetGML,'r'):
		if l.find('node')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
			if idT in sIdsIntrst:
				hldr.append('  node [\n    id %s\n'%idT)
		elif hldr:
			hldr.append(l)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest 
	hldr = []
	sIdTInEdgs = set()
	cntFnlEdgs = 0
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst and trgt in sIdsIntrst and \
				{sorted([src,trgt])}.issubset(sShrdVldEdgs):
				cntFnlEdgs+=1
				sIdTInEdgs.update({src,trgt})
				frstNtwrk = '1.0'
			else:
				frstNtwrk = '0.0'
			hldr.append \
			('  edge [\n    source %s\n    target %s\n    frstNtwrk %s\n'% \
			(src,trgt,frstNtwrk))
		elif hldr:
			hldr.append(l)
	if hldr:#last case
		G_targetGMLOut.write(''.join(hldr))
	G_targetGMLOut.write(']')
	G_targetGMLOut.close()
	#----------------------------
	# report final results
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges and %s nodes of interest are'% \
		(cntFnlEdgs,len(sIdTInEdgs)), \
		'were overlaping and retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Make a matrix from a GML graph file with a graph of interest
########################################################
def readGMLtMtrx(G_targetGML,aVlsANms=None,aVlsBNms=None,vrbse=True):
	"""
	Input: G_targetGML is the input graph of interest (from which a 
	subset is going to be obtained) in '.gml' format. Optionally, 
	aVlsANms is an array of names for rows. aVlsBNms is an array of 
	names for columns. If vrbse is True, the log for the run is going to 
	be printed.
	Output: aVlsAVlsB is a 2x2 matrix of values of size len(aVlsANms) x
	len(aVlsBNms). Values for rows must be in the same order as names in
	aVlsANms, and values for columns in the same order as names in 
	aVlsBNms. Optionally, aVlsAVlsB and aVlsANms are the array of names.
	ovrAllogRun is the log of the run.
	NOTE: If aVlsANms and aVlsBNms are not provided the method will 
	make a new array from the sorted names in the input file.
	"""
	#----------------------------
	#~ get gml information for nodes of interest and make intial matrix
	dIdsNm = {}
	dNmIds = {}
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			assert not dNmIds.has_key(nodeNm)
			dIdsNm[idT] = nodeNm
			dNmIds[nodeNm] = idT
	srtdNms = sorted(set(dIdsNm.values()))
	if aVlsANms is not None and aVlsBNms is not None:
		dNmsPosinaVlsANms = dict([((nmA,pos) for pos,nmA in \
		enumerate(aVlsANms))])
		dNmsPosinaVlsBNms = dict([((nmB,pos) for pos,nmB in \
		enumerate(aVlsBNms))])
		sVlsANms,sVlsBNms = set(aVlsANms),set(aVlsBNms)
		sNdsInGrph = set(srtdNms)
		assert sNdsInGrph.issubset(sVlsANms) and sNdsInGrph.issubset \
		(sVlsBNms)
		lenVlsANms,lenVlsBNms = len(aVlsANms),len(aVlsBNms)
		aVlsAVlsB = zeros((lenVlsANms,lenVlsBNms),dtype=float32)
		tstBthlst = True
	else:
		lenVlsNms = len(srtdNms)
		dNmsPosinaVlsNms = dict([((nm,pos) for pos,nm in \
		enumerate(srtdNms))])
		aVlsAVlsB = zeros((lenVlsNms,lenVlsNms),dtype=float32)
		tstBthlst = False
	aVlsAVlsB.fill(nan)
	#----------------------------
	# report preliminary
	ovrAllogRun = []
	try:
		if tstBthlst:
			raise exceptions.CelleryWarningObjct \
			('%s nodes for A values and %s nodes for B values'% \
			(lenVlsANms,lenVlsBNms),'were found in the input file')
		else:
			raise exceptions.CelleryWarningObjct \
			('%s nodes'%lenVlsNms,'were found in the input file')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ get gml information for edges of interest and fill output matrix
	hldr = []
	fnlEdgs,setFnlNds = 0,set()
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
		elif l.find('weight')>-1:
			wght = float32(l.split('weight ')[1].splitlines()[0])
			if tstBthlst:
				if src in sVlsANms and trgt in sVlsBNms:
					posValA = dNmsPosinaVlsANms[dIdsNm[src]]
					posValB = dNmsPosinaVlsBNms[dIdsNm[trgt]]
					aVlsAVlsB[posValA][posValB] = wght
					fnlEdgs+=1
					setFnlNds.update({posValA,posValB})
				elif src in sVlsBNms and trgt in sVlsANms:
					posValA = dNmsPosinaVlsANms[dIdsNm[trgt]]
					posValB = dNmsPosinaVlsBNms[dIdsNm[src]]
					aVlsAVlsB[posValA][posValB] = wght
					fnlEdgs+=1
					setFnlNds.update({posValA,posValB})
				else:
					raise exceptions.CelleryExceptionObjct \
					('%s nodes for A values and %s nodes for B values')
			else:
				posValA = dNmsPosinaVlsNms[dIdsNm[trgt]]
				posValB = dNmsPosinaVlsNms[dIdsNm[src]]
				aVlsAVlsB[posValA][posValB] = wght
				fnlEdgs+=1
				setFnlNds.update({posValA,posValB})
	#----------------------------
	# final report
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes and %s edges'%(len(setFnlNds),fnlEdgs), \
		'were retrieved from the input file')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	if tstBthlst:
		return aVlsAVlsB,aVlsANms,aVlsBNms,ovrAllogRun
	else:
		aVlsANms = array(srtdNms)
		aVlsBNms = array(srtdNms)
		return aVlsAVlsB,aVlsANms,aVlsBNms,ovrAllogRun


########################################################
#~ Return a subgraph made up by a subset of nodes of interest using GML 
# graph files
########################################################
def rtrnGMLSubGrph(sNds,G_targetGML,G_targetGMLOut,vrbse=True):
	"""
	Input: sNds is a set of nodes of interest to substract from a graph. 
	G_targetGML is the input graph of interest (from which a subset is 
	going to be obtained) in '.gml' format. G_targetGMLOut is the output 
	file to write the subset graph in '.gml' format. If vrbse is True, 
	the log forthe run is going to be printed.
	Output: G_targetGMLOut is the output file with a subgraph of nodes 
	of interest and the edges between them. ovrAllogRun is the log of 
	the run.
	NOTE: Assumes undirected graphs for input and output.
	"""
	#----------------------------
	#~ get nodes of interest
	sIdsIntrst = set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			if nodeNm in sNds:
				sIdsIntrst.add(idT)
	#----------------------------
	#~ set parameters
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest'%len(sIdsCnntng), \
		'are going to be retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ get gml information for nodes of interest 
	G_targetGMLOut = open(G_targetGMLOut,'w')
	G_targetGMLOut.write('graph [\n')
	hldr = []
	cntFnlNds = 0
	for l in open(G_targetGML,'r'):
		if l.find('node')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
			if idT in sIdsIntrst:
				cntFnlNds+=1
				hldr.append('  node [\n    id %s\n'%idT)
		elif hldr:
			hldr.append(l)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest 
	hldr = []
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst and trgt in sIdsIntrst:
				hldr.append('  edge [\n    source %s\n    target %s\n'% \
				(src,trgt))
		elif hldr:
			hldr.append(l)
	if hldr:#last case
		G_targetGMLOut.write(''.join(hldr))
	G_targetGMLOut.write(']')
	G_targetGMLOut.close()
	#----------------------------
	#~ report final values
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest'%cntFnlNds,'were retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Return a subgraph made up by a subset of nodes of interest and their
# immediate neighbors using GML graph files
########################################################
def rtrnGMLSubGrphClstNghbrs(sNds,G_targetGML,G_targetGMLOut,vrbse=True):
	"""
	Input: sNds is a set of nodes of interest to substract from a graph. 
	G_targetGML is the input graph of interest (from which a subset is 
	going to be obtained) in '.gml' format. G_targetGMLOut is the output 
	file to write the subset graph in '.gml' format. If vrbse is True, 
	the log forthe run is going to be printed.
	Output: G_targetGMLOut is the output file with a subgraph of nodes 
	of interest, their immediate neighbors, and the edges between them.
	ovrAllogRun is the log of the run.
	NOTE: Assumes undirected graphs for input and output.
	"""
	#----------------------------
	#~ get nodes of interest, and node with edges shared with the nodes 
	# of interest
	sIdsIntrst,sIdsCnntng = set(),set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			if nodeNm in sNds:
				sIdsIntrst.add(idT)
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst:
				sIdsCnntng.add(trgt)
			elif trgt in sIdsIntrst:
				sIdsCnntng.add(src)
	#----------------------------
	#~ set parameters
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest and %s immediate neighbors'% \
		(len(sIdsIntrst),len(sIdsCnntng)),'are going to be retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	sIdsIntrst.update(sIdsCnntng)
	#----------------------------
	#~ get gml information for nodes of interest 
	G_targetGMLOut = open(G_targetGMLOut,'w')
	G_targetGMLOut.write('graph [\n')
	hldr = []
	for l in open(G_targetGML,'r'):
		if l.find('node')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
			if idT in sIdsIntrst:
				hldr.append('  node [\n    id %s\n'%idT)
		elif hldr:
			hldr.append(l)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest 
	hldr = []
	sFnlNds = set()
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst and trgt in sIdsIntrst:
				hldr.append('  edge [\n    source %s\n    target %s\n'% \
				(src,trgt))
				sFnlNds.update({src,trgt})
		elif hldr:
			hldr.append(l)
	if hldr:#last case
		G_targetGMLOut.write(''.join(hldr))
	G_targetGMLOut.write(']')
	G_targetGMLOut.close()
	#----------------------------
	#~ report final values
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest and immediate neighbors'%len(sFnlNds), \
		'were retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Return a subgraph made up by a subset of edges of interest using GML 
# graph files
########################################################
def rtrnGMLSubGrphEdgs(sTplSrtdEdgsFrstNtwrk,G_targetGML,G_targetGMLOut, \
	vrbse=True):
	"""
	Input: sTplSrtdEdgsFrstNtwrk is a set of sorted pair-nodes 
	representing edges of interest to substract from a graph. 
	G_targetGML is the input graph of interest (from which a subset is 
	going to be obtained) in '.gml' format. G_targetGMLOut is the output 
	file to write the subset graph in '.gml' format. If vrbse is True, 
	the log forthe run is going to be printed.
	Output: G_targetGMLOut is the output file with a subgraph of edges 
	of interest. ovrAllogRun is the log of the run.
	NOTE: Assumes undirected graphs for input and output.
	"""
	#----------------------------
	#~ get nodes of in the graph
	sNds = set()
	for ndA,ndB in sTplSrtdEdgsFrstNtwrk:
		sNds.update({ndA,ndB})
	sIdsIntrst = set()
	dNodeNmIdT = {}
	sNodeNms = set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			if nodeNm in sNds:
				sIdsIntrst.add(idT)
				dNodeNmIdT[nodeNm] = idT
				sNodeNms.add(nodeNm)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get nodes of interest
	sShrdVldEdgs = []
	while len(sTplSrtdEdgsFrstNtwrk)>0:
		ndA,ndB = sTplSrtdEdgsFrstNtwrk.pop()
		if ndA in sNodeNms and ndB in sNodeNms:
			idTA,idTB = dNodeNmIdT[ndA],dNodeNmIdT[ndB]
			sShrdVldEdgs.append(sorted([idTA,idTB]))
	#----------------------------
	# report preliminary
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges and %s nodes of interest are'% \
		(len(sTplSrtdEdgsFrstNtwrk),len(sNodeNms)), \
		'going to be searched')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ get gml information for nodes of interest
	G_targetGMLOut = open(G_targetGMLOut,'w')
	G_targetGMLOut.write('graph [\n')
	hldr = []
	for l in open(G_targetGML,'r'):
		if l.find('node')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
			if idT in sIdsIntrst:
				hldr.append('  node [\n    id %s\n'%idT)
		elif hldr:
			hldr.append(l)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest 
	hldr = []
	sIdTInEdgs = set()
	cntFnlEdgs = 0
	for l in open(G_targetGML,'r'):
		if l.find('edge')>-1:
			if hldr:
				G_targetGMLOut.write(''.join(hldr))
			hldr = []
		elif l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
			if src in sIdsIntrst and trgt in sIdsIntrst and \
				{sorted([src,trgt])}.issubset(sShrdVldEdgs):
				cntFnlEdgs+=1
				sIdTInEdgs.update({src,trgt})
				hldr.append('  edge [\n    source %s\n    target %s\n'% \
				(src,trgt))
		elif hldr:
			hldr.append(l)
	if hldr:#last case
		G_targetGMLOut.write(''.join(hldr))
	G_targetGMLOut.write(']')
	G_targetGMLOut.close()
	#----------------------------
	# report final results
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges and %s nodes of interest are'% \
		(cntFnlEdgs,len(sIdTInEdgs)), \
		'were overlaping and retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Return a set of nodes a list of weighted edges from a GML graph file
########################################################
def readGMLtSNdslWghtdEdgs(G_targetGML,vrbse=True):
	"""
	Input: G_targetGML is the input graph of interest (from which a 
	subset is going to be obtained) in '.gml' format. If vrbse is True, 
	the log for the run is going to be printed.
	Output: sNodeNms is the set of nodes in the input graph file. 
	lWghtEdgs is the list of edges with weight. ovrAllogRun is the log 
	of the run.
	"""
	#----------------------------
	#~ get gml information for nodes of interest
	G_targetGMLOut = open(G_targetGMLOut,'w')
	sIdsIntrst = set()
	dIdTNodeNm = {}
	sNodeNms = set()
	for l in open(G_targetGML,'r'):
		if l.find('id ')>-1:
			idT = l.split('id ')[1].splitlines()[0]
		elif l.find('label')>-1:
			nodeNm = l.split('"')[1]
			sIdsIntrst.add(idT)
			dIdTNodeNm[idT] = nodeNm
			assert nodeNm not in sNodeNms#not repeated ids in each node
			sNodeNms.add(nodeNm)
		elif l.find('edge')>-1:
			break
	#----------------------------
	#~ get gml information for edges of interest and fill output matrix
	lWghtEdgs = []
	for l in open(G_targetGML,'r'):
		if l.find('source')>-1:
			src = l.split('source ')[1].splitlines()[0]
		elif l.find('target')>-1:
			trgt = l.split('target ')[1].splitlines()[0]
		elif l.find('weight')>-1:
			wght = float32(l.split('weight ')[1].splitlines()[0])
			lWghtEdgs.append([wght,dIdTNodeNm[src],dIdTNodeNm[trgt]])
	#----------------------------
	# final report
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes and %s edges'%(len(sNodeNms),lWghtEdgs), \
		'were retrieved from the input file')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return sNodeNms,lWghtEdgs,ovrAllogRun


#----------------------------
#~ networkX based methods
#----------------------------

########################################################
#~ Annotate edges with information of interest in a networkX graph
########################################################
def anntNtwrkXEdgs(G_target,dEdgesAnnt,anntNm,vrbse=True):
	"""
	Input: G_target is the graph of interest in networkx. dEdgesAnnt is 
	the dictionary with edges as keys and probabilities/annotations as 
	values. If vrbse is True, the log forthe run is going to be printed.
	Output: G_target with edges in the keys dEdgesAnnt annotated. If the
	edges are not present in the keys of dEdgesAnnt, the are anotated as
	'N'. ovrAllogRun is the log of the run.
	"""
	#----------------------------
	# Method to parse the input annotations
	ovrAllogRun = []#log holder
	sEdgs = set(dEdgesAnnt.keys())
	lensEdgs = len(sEdgs)
	nonAnntd = 0
	for nodA,nodB in sEdgs:
		annt = dEdgesAnnt.get((nodA,nodB),False)
		if not annt:
			nonAnntd+=1
			annt = 'N'
		G_target[nodA][nodB][anntNm] = annt
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges are present in the input graph'%lensEdgs, \
		'%s were annotated'%lensEdgs-nonAnntd)
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Annotate edges in a target networkX graph in absent in a reference 
# one
########################################################
def anntNtwrkXEdgsNotInRef(G_target,G_ref,vrbse=True):
	"""
	Input: G_target is the graph of interest, G_ref is the graph of 
	reference. If vrbse is True, the log forthe run is going to be 
	printed.
	Output: G_target (is a graph in networkx format) with edges not 
	present in G_ref colored in red. ovrAllogRun is the log of the run.
	"""
	#----------------------------
	# remove edges
	sMrgEdgs = set(G_target.edges()).union(set(G_ref.edges()))
	crldBlckRef,crldRedRef = 0,0
	crldBlckTrgt,crldRedTrgt = 0,0
	for nodA,nodB in G_ref.edges():
		if  G_target.has_edge(nodA,nodB):
			G_target[nodA][nodB]['color'] = 'black'
			crldBlckTrgt+=1
		else:
			G_target.add_edge(nodA,nodB,{'weight':G_ref[nodA][nodB] \
			['weight'],'color':'red'})
			crldRedTrgt+=1
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('and %s edges were colored black'%crldBlckTrgt, \
		'and %s red for the target graph'%crldRedTrgt)
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Annotate edges in a target networkX graph in absent in a reference 
# one
########################################################
def anntNtwrkXEdgsNotInRefNVicvrs(G_target,G_ref,vrbse=True):
	"""
	Input: G_target is the graph of interest, G_ref is the graph of 
	reference. If vrbse is True, the log forthe run is going to be 
	printed.
	Output: G_target (is a graph in networkx format) with edges not 
	present in G_ref colored in red and viceversa. ovrAllogRun is the 
	log of the run.
	"""
	#----------------------------
	# remove edges
	sMrgEdgs = set(G_target.edges()).union(set(G_ref.edges()))
	crldBlckRef,crldRedRef = 0,0
	crldBlckTrgt,crldRedTrgt = 0,0
	for nodA,nodB in sMrgEdgs:
		if  G_ref.has_edge(nodA,nodB):
			G_ref[nodA][nodB]['color'] = 'black'
			crldBlckRef+=1
		else:
			G_ref.add_edge(nodA,nodB,{'weight':G_target[nodA][nodB] \
			['weight'],'color':'red'})
			crldRedRef+=1
		if  G_target.has_edge(nodA,nodB):
			G_target[nodA][nodB]['color'] = 'black'
			crldBlckTrgt+=1
		else:
			G_target.add_edge(nodA,nodB,{'weight':G_ref[nodA][nodB] \
			['weight'],'color':'red'})
			crldRedTrgt+=1
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s edges were colored black and %s red for the reference graph' \
		%(crldBlckRef,crldRedRef), \
		'and %s edges were colored black and %s red for the target graph' \
		%(crldBlckTrgt,crldRedTrgt))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun


########################################################
#~ Make a networkX graph from a gml file
########################################################
def mkNtwrkXFrmGML(G_targetGML,aVlsANms=None,aVlsBNms=None,vrbse=True):
	"""
	Input: G_targetGML is the input graph of interest (from which a 
	subset is going to be obtained) in '.gml' format. Optionally, 
	aVlsANms is an array of names for rows. aVlsBNms is an array of 
	names for columns. If vrbse is True, the log for the run is going to 
	be printed.
	Output: A graph in networkX format with cell values as edge weights 
	(annotated under the label "weight"). ovrAllogRun is the log of the 
	run.
	NOTE: If aVlsANms and aVlsBNms are not provided the method will 
	make a new array from the sorted names in the input file.	
	"""
	#----------------------------
	#~ obtain a matrix from a gml file
	aVlsAVlsB,aVlsANms,aVlsBNms,ovrAllogRun = readGMLtMtrx(G_targetGML, \
	aVlsANms,aVlsBNms,vrbse)
	#----------------------------
	#~ set initial parameters
	G = nx.Graph()
	lenVlsA,lenVlsB = aVlsAVlsB.shape
	assert len(aVlsANms)==lenVlsA and len(aVlsBNms)==lenVlsB
	#----------------------------
	#~ set initial parameters
	cntEdgs = 0
	sNdsFnl = set()
	for vlsBPos in xrange(lenVlsB):
		vlsBNm = aVlsBNms[vlsBPos]
		for vlsAPos in xrange(lenVlsA):
			vlsANm = aVlsANms[vlsAPos]
			if undrctd and vlsBNm==vlsANm:
				pass
			else:
				prwsVal = aVlsAVlsB[vlsAPos,vlsBPos]
				if not isnan(dtPrws):
					cntEdgs+=1
					sNdsFnl.update({vlsANm,vlsBNm})
					G.add_edge(vlsANm,vlsBNm,{'weight':str(prwsVal)})
	#----------------------------
	# report final results
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges was build'%(cntEdgs,'for %s nodes.'% \
		len(sNdsFnl)))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G,ovrAllogRun


########################################################
#~ Make a networkX graph from a matrix with cells representing edges
########################################################
def mkNtwrkXFrmMtrx(aVlsAVlsB,aVlsANms,aVlsBNms,undrctd=False,vrbse=True):
	"""
	Input: aVlsAVlsB is a 2x2 matrix of values of size len(aVlsANms) x
	len(aVlsBNms). Values for rows must be in the same order as names in
	aVlsANms, and values for columns in the same order as names in 
	aVlsBNms. Optionally, if undrctd is True only values for the first 
	(sorted) names is going to be considered. If vrbse is True, the log 
	for the run is going to be printed.
	Output: G is a graph in networkX format with cell values as edge 
	weights (annotated under the label "weight"). ovrAllogRun is the log 
	of the run.
	"""
	#----------------------------
	#~ set initial parameters
	G = nx.Graph()
	lenVlsA,lenVlsB = aVlsAVlsB.shape
	assert len(aVlsANms)==lenVlsA and len(aVlsBNms)==lenVlsB
	#----------------------------
	#~ set initial parameters
	cntEdgs = 0
	sNdsFnl = set()
	for vlsBPos in xrange(lenVlsB):
		vlsBNm = aVlsBNms[vlsBPos]
		for vlsAPos in xrange(lenVlsA):
			vlsANm = aVlsANms[vlsAPos]
			if undrctd and vlsBNm==vlsANm:
				pass
			else:
				prwsVal = aVlsAVlsB[vlsAPos,vlsBPos]
				if not isnan(dtPrws):
					cntEdgs+=1
					sNdsFnl.update({vlsANm,vlsBNm})
					G.add_edge(vlsANm,vlsBNm,{'weight':str(prwsVal)})
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges was build'%(cntEdgs,'for %s nodes.'% \
		len(sNdsFnl)))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G,ovrAllogRun


########################################################
#~ Make a networkX graph from a matrix with cells representing edges
########################################################
def mkNtwrkXFrmTSV(inTSV,dlmtr='\t',vrbse=True):
	"""
	Input: inTSV is a file in txt format with values in cells between 
	different rows (values A) and columns (values B). dlmtr is the 
	delimiter to split the file. If vrbse is True, the log for the run 
	is going to be printed. 
	Output: G_target is a graph in networkX format with cell values as 
	edge weights (annotated under the label "weight"). aVlsANms is an 
	array of names for rows. aVlsBNms is an array of names for columns. 
	ovrAllogRun is the log of the run.
	NOTE: inTSV must have a header with names of variables (i.e. genes/
	isoforms) in columns and an additional column (in position 0) with 
	names of variables (i.e. genes/isoforms) in rows.
	"""
	#----------------------------
	# set initial variables
	assert prctWght or trhldEdgs
	G_target = nx.Graph()
	intbl = np.loadtxt(inTSV,delimiter=dlmtr,dtype=str)
	aVlsANms = []
	hdr = True
	nEdges = 0
	for el in intbl:
		if hdr:
			nVlsB = len(el)
			aVlsBNms = array(el[1:])
			dPosVlB=dict([(pos,VlB) for pos,VlB in enumerate(el)])
			hdr = False
		else:
			nonNm = False
			for echPos in xrange(nVlsB):
				if nonNm:
					val = el[echPos]
					vlBNm = dPosVlB[echPos]
					if val not in {'nan',''}:
						G_target.add_edge(vlANm,vlBNm,{'weight': \
						float32(val)})
						nEdges+=1
				else:
					vlANm = el[echPos]
					aVlsANms.append(vlANm)
					nonNm = True
	aVlsANms = array(aVlsANms)
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges was retrieved'%(nEdges,
		'between %s nodes.'% len(aVlsANms)+len(aVlsBNms)))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G_target,aVlsANms,aVlsBNms,ovrAllogRun


########################################################
#~ Remove edges in a target networkX graph in absent in a reference one
########################################################
def rmvNtwrkXEdgsByWght(G_target,prctCut=False,trhldCut=False,vrbse=True):
	"""
	Input: G_target is the graph of interest. prctCut is a percentile 
	(0-100) under which edges are going to be cut. trhldCut is a weight 
	value under which edges are going to be cut. If vrbse is True, the 
	log for the run is going to be printed. 
	Output: G_targetCut (is a graph in networkx format) representing the
	input graph without the edges filtered out by trhldCut or vrbse 
	parameters. aVlsANms is an array of names for rows. aVlsBNms is an 
	array of names for columns. ovrAllogRun is the log of the run.
	"""
	#----------------------------
	# test for input parameters
	assert (trhldCut and not prctCut) or (prctCut and not trhldCut)
	G_targetCut = G_target.copy()
	#----------------------------
	# preliminary report
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges'%G_targetCut. \
		number_of_edges(),'and %s nodes will be cutted'%G_targetCut. \
		number_of_nodes())
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	# select edges following the input parameters
	if trhldEdgs:
		edges_del = [(edgA,edgB) for edgA,edgB in G_targetCut.edges() if \
		G_targetCut[edgA][edgB]['weight']>trhldEdgs]#like in a probability
	else:
		trhldEdgs = percentile([edge[2]['weight'] for edge in \
		G_targetCut.edges(data=True)],prctWght)
		edges_del = [(edgA,edgB) for edgA,edgB in G_targetCut.edges() if \
		G_targetCut[edgA][edgB]['weight']<trhldEdgs]#like in a score
	G_targetCut.remove_edges_from(edges_del)
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges'%G_targetCut. \
		number_of_edges(),'and %s nodes were retrieved'%G_targetCut. \
		number_of_nodes())
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G_targetCut,aVlsANms,aVlsBNms,ovrAllogRun


########################################################
#~ Remove edges in a target networkX graph in absent in a reference one
########################################################
def rmvNtwrkXEdgsNotInRef(G_target,G_ref,vrbse=True):
	"""
	Input: G_target is the graph of interest. G_ref is the graph of 
	reference.
	Output: G_targetCut (is a graph in networkx format) representing the
	input graph without the edges not present in G_ref. ovrAllogRun is 
	the log of the run.
	"""
	#----------------------------
	# preliminary report
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges'%G_targetCut. \
		number_of_edges(),'and %s nodes will be cutted'%G_targetCut. \
		number_of_nodes())
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	# remove edges
	G_targetCut = G_target.copy()
	for nodA,nodB in G_targetCut.edges():
		if not G_ref.has_edge(nodA,nodB):
			G_targetCut.remove_edge(nodA,nodB)
	for nod in G_targetCut.nodes():
		if G_targetCut.degree(nod)==0:
			G_targetCut.remove_node(nod)
	#----------------------------
	# report final results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A graph with %s edges'%G_targetCut. \
		number_of_edges(),'and %s nodes were retrieved'%G_targetCut. \
		number_of_nodes())
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G_targetCut,ovrAllogRun


########################################################
#~ Write networkX graph to a GML graph file
########################################################
def writeNtwrkXToGML(G_target,G_targetGMLOut):
	"""
	Input: G_target is the graph of interest. G_targetGMLOut is the 
	output file to write the subset graph in '.gml' format.
	Output: G_targetGMLOut is the output file to write the subset graph 
	in '.gml' format.
	"""
	nx.write_gml(G_target,G_targetGMLOut)
	#~ return 0


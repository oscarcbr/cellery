#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  __classes__.py part of cellery (ceRNAs linking inference)
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
Define gene class with pointers to different propierties of the genes
"""


########################################################
#~ Import libraries.
########################################################
from cellery import exceptions
from numpy import array,zeros,float32


########################################################
#~ Define gene class
########################################################
class gene(object):
	"""
	Object to represent gene/lncRNA. It has pointers 1) "name" for the 
	gene/lncRNA name, 2) "len" for the length of gene/lncRNA size (in 
	nt), 3) "pos" for the position of the gene/lncRNA in the list of 
	gene objects (the sorted names of all genes of iterest), 4) 
	"cmmNm" common/external gene name, 5) "strnd" is the strand of the 
	gene, 6) aIntrvls is the array of intervals for exons/sequences of
	each gene. 7) "aMirnaNms" for an array of the sorted names of miRNAs 
	(headers) of interest. 8) "aMirnaCnts" that is an array of MRE 
	counts for a gene interest predicted by TargetScan. Further pointers 
	to other program results are included as follows: aMirnaScrs for 
	context+ scores, aMirnaRNAhCnts for RNAhybrid counts, aMirnaRNAhEngy
	for RNAhybrid (minimal) free energy, aMirnaMirndCnts for miRanda 
	counts, aMirnaMirndScrs for miRanda scores, aMirnaMirndEngy for 
	miRNA (minimal) free energy, aMirnaSVMicroCnts for SVMicro counts, 
	aMirnaSVMicroScrs for SVMicro scores, aMirnaTrgtMnrCnts for 
	TargetMiner, aMirnaPITACnts for PITA counts, aMirnaPITACnts for 
	PITA scores, aMirnaMirMapScrs for mirMap scores, aMirnaMirMapCnts 
	for mirMap counts, aMirnaWalk for miRWalk prediction status.
	"""
	__slots__ = ('name','len','pos','chr','cmmNm','strnd','aIntrvls', \
	'aMirnaNms','aMirnaCnts','aMirnaScrs','aMirnaRNAhCnts', \
	'aMirnaRNAhEngy','aMirnaMirndCnts','aMirnaMirndScrs', \
	'aMirnaMirndEngy','aMirnaSVMicroCnts','aMirnaSVMicroScrs', \
	'aMirnaMirndEngy','aMirnaSVMicroCnts','aMirnaSVMicroScrs', \
	'aMirnaTrgtMnrCnts','aMirnaPITACnts','aMirnaPITAScrs', \
	'aMirnaMirMapScrs','aMirnaMirMapCnts','aMirnaWalk')
	def __init__(self,name,srtdMirnaNms):
		self.name = name
		self.len = False
		self.pos = False
		self.chr = set()
		self.cmmNm = False
		self.strnd = False
		self.aIntrvls = False
		cntMirnas = len(srtdMirnaNms)
		self.aMirnaNms = array(srtdMirnaNms)
		self.aMirnaCnts = zeros(cntMirnas,dtype=float32)
		#~ self.aMirnaScrs = zeros(cntMirnas,dtype=float32)#void contxt+
		self.aMirnaRNAhCnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaRNAhEngy = zeros(cntMirnas,dtype=float32)
		self.aMirnaMirndCnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaMirndScrs = zeros(cntMirnas,dtype=float32)
		self.aMirnaMirndEngy = zeros(cntMirnas,dtype=float32)
		self.aMirnaSVMicroCnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaSVMicroScrs = zeros(cntMirnas,dtype=float32)
		self.aMirnaTrgtMnrCnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaPITACnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaPITAScrs = zeros(cntMirnas,dtype=float32)
		self.aMirnaMirMapScrs = zeros(cntMirnas,dtype=float32)
		self.aMirnaMirMapCnts = zeros(cntMirnas,dtype=float32)
		self.aMirnaWalk = zeros(cntMirnas,dtype=bool)


########################################################
#~ Define node class
########################################################
class node(object):
	"""
	Object to represent nodes. It has pointers 1) "name" for the node 
	name, 2) "w" for the weight of node size (is a float32 value),
	"""
	__slots__ = ('w', 'name')
	def __init__(self,name):
		self.name = name
		self.w = float32(0)


########################################################
#~ Define disjoint-set class
########################################################
class disjntSet(object):
	"""
	Object to represent disjoint-set data structure as defined by Cormen 
	et al. 2013 page 561. It keeps track of a set of elements 
	partitioned into a number of disjoint (nonoverlapping) subsets.
	"""
	def __init__(self,lX):
		self.clltn = []
		for x in lX:
			self.clltn.append({x})
	#----------------------------
	def makeSet(self,x):
		"""
		Method to make a new set disjoint set. Creates the pointer 
		clltn.
		"""
		prsnt = False
		if self.clltn:
			for eSet in self.clltn:
				if x in eSet:
					prsnt = True
					break
		if not prsnt:
			self.clltn.append({x})#if not present create a new set
	#----------------------------
	def djntUnion(self,x,y):
		"""
		Unites the dynamic sets that contain x and y, say Sx and Sy, 
		into a new set that is the union of these two sets.
		"""
		jnd = False
		xSet = False#to assert x is present
		ySet = False#to assert y is present
		newClltn = []
		if self.clltn:
			while len(self.clltn)>0:
				cSet = self.clltn.pop()
				appnd = True
				if x in cSet:
					xcSet = cSet
					xSet = True
					appnd = False
				if y in cSet:
					ycSet = cSet
					ySet = True
					appnd = False
				if appnd:
					newClltn.append(cSet)
				elif xSet and ySet:
					newClltn.append(xcSet.union(ycSet))
					newClltn.extend(self.clltn)
					break
		if not xSet and not ySet:
			raise exceptions.CelleryExceptionObjct \
			('%s or %s have not been defined'%(x,y))
		else:
			self.clltn = newClltn
	#----------------------------
	def findSet(self,x):
		"""
		Returns a pointer to the representative of the (unique) set 
		containing x.
		"""
		pos = -1
		for eSet in self.clltn:
			pos+=1
			if x in eSet:
				break
		return pos

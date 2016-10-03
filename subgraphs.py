#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  subgraphs.py part of cellery (ceRNAs linking inference)
#  
#  Copyright 2015 Oscar Bedoya Reina <obedoya@igmm-linux-005>
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
Libraries to randomize networks and subgraphs.
"""

########################################################
#~ Import libraries.
########################################################
#----------------------------
#Import matplotlib libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#----------------------------
#Import other libraries
from cellery import __class__
from cellery import exceptions
from cellery import networks

from itertools import combinations,count
from multiprocessing import Queue,Process
from numpy import array,float32,float64,empty,exp,ma,nan,percentile
from random import choice
from scipy.stats import bernoulli
from string import upper

import argparse
import networkx as nx
import numpy as np
import os


########################################################
#~ Return important nodes (i.e. nodes that have a high weighted degree)
# sharing a similar structure (neighbors of neighbors) by two networks.
########################################################
def slctImprtntNds(G_ntwrkxA,G_ntwrkxB,rho=0.7,theta=0.7,nThrd=5, \
	vrbse=True):
	"""
	Input: G_ntwrkxA is one graph of interest, for instance generated 
	from the endogenous competition data. G_ntwrkxB is one graph of 
	interest, for instance generated from the expression data. rho is 
	the degree of approximation orthe percentage of neighbors of a query 
	node that can have no corresponding matches in the neighborhood of a 
	database node, and theta is equivalent to rho, just that it is 
	applicable to the sum of the weighted edges between neighbors. 
	Neighbor connections are the number of edges between the neighbors. 
	If vrbse is True, the log forthe run is going to be printed.
	Output: A list of nodes that satifies the conditions of rho and 
	theta, and pointers to the quality and names of all nodes.
	Note: This module runs in multithread, the parameters nThrd is the
	number of cores to be used by the method.
	"""
	#----------------------------
	#~ Core Method
	def core(G_ntwrkxA,G_ntwrkxB,nod,rho,theta):
		nodQDgreeGed,nbMissGed,nbcMissGed,nbConnctinGed=clctScrs \
		(G_ntwrkxA,nod,rho,theta)
		#G_ntwrkxA==query graph
		nodQDgreeGce,nbMissGce,nbcMissGce,nbConnctinGce=clctScrs \
		(G_ntwrkxB,nod,rho,theta)		
		nbMissActl,nbcMissActl=calcMiss(G_ntwrkxA,G_ntwrkxB,nod)
		if nbMissActl<=nbMissGed and nodQDgreeGed>0:#if True when 
			#nbMissGed==0, it implies that all the neighbors of "node" 
			# in G_ntwrkxA are present in G_ntwrkxB.
			fNb=nbMissActl/float(nodQDgreeGed)
			if nbConnctinGed>0:
				fNbc=nbcMissActl/float(nbConnctinGed)
			else:
				fNbc=0.0
			#Therefore, we amortize fNbc by the number of missing 
			#neighbors nbmiss. The value of w falls between 0 and 2. We 
			#substract this value from 2, so that higher w value means a 
			#better node match.		
			if nbMissActl==0:
				w=2-fNbc#two is the maximum weight between nodes
			else:
				w=2-(fNb+(fNbc/nbMissActl))
			nodObj=__class__.node(nod)
			nodObj.w=w
			return nodObj
		else:
			return False
	#----------------------------
	#~ Multithread method
	def mltthrd(qInJobs,qOutRslts,G_ntwrkxA,G_ntwrkxB,rho,theta):
		for nod in iter(qInJobs.get,'STOP'):
			nodObj=core(G_ntwrkxA,G_ntwrkxB,nod,rho,theta)
			qOutRslts.put(nodObj)
	#----------------------------
	#~ Make a set of intersecting nodes
	Mimp=[]
	sShrdNods=set(G_ntwrkxA.nodes()).intersection(set \
	(G_ntwrkxB.nodes()))
	nShrdNods=len(sShrdNods)
	#----------------------------
	# report preliminary
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('%s nodes of interest are'%nShrdNods,'going to be searched')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	#~ Run
	qInJobs = Queue()
	qOutRslts = Queue()
	for nod in sShrdNods:
		qInJobs.put(nod)
	for p in range(nThrd):
		Process(target=mltthrd,args=(qInJobs,qOutRslts,G_ntwrkxA, \
		G_ntwrkxB,rho,theta)).start()
	c=0
	for p in range(nShrdNods):
		c+=1
		nodObj=qOutRslts.get()
		if nodObj:
			Mimp.append(nodObj)
		if vrbse and c%100==0:
			print 'Running calculations on node %s out of %s'% \
			(c,nShrdNods)
	for p in range(nThrd):
			qInJobs.put('STOP')
	#----------------------------
	# final report
	try:
		raise exceptions.CelleryWarningObjct \
		('%s important shared nodes'%len(Mimp),'were retrieved')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return Mimp,ovrAllogRun


########################################################
#~ Compute the maximum number of neighbor missed nodes following rho and 
# theta. These values are required to intersect networks. 
########################################################
def clctScrs(G,nod,rho,theta):
	"""
	Input: G is a grpah of interest, and nod is a node of interest, rho
	is the degree of approximation or the percentage of neighbors of a 
	query node that can have no corresponding matches in the 
	neighborhood of a database node, and theta is equivalent to rho but 
	for the connections of neighbors.
	Output: nbMiss is the number of missing neighbors in the node match, 
	and nbcMiss is the number of missing connections, nodeQDgree is the
	degree of the node of interest, nbConnctin is the sum of the 
	weighted edges bewteen the neighbors.
	"""
	nodQDgree = G.degree(nod,'weight')#weighted node degree
	nbMiss = rho*nodQDgree#nodQDgree neighbors of the query node
	#can be missing in the match to a database node
	nbConnctin = 0.0
	lNodNeigbr = G[nod].keys()
	for nodNeigbrA,nodNeigbrB in combinations(lNodNeigbr,2):
		if G.has_edge(nodNeigbrA,nodNeigbrB):
			nbConnctin+=G[nodNeigbrA][nodNeigbrB]['weight']
	nbcMiss = nbConnctin*theta#nbcMiss nodes redifined from the original 
	#paper: (nbMiss*(nbMiss-1))/(2+(nbMiss*(nodQDgree-nbMiss))) "all 
	#connect to each other, and also connect to all of the remaining 
	#(nodQDgree-nbMiss) nodes. This must be changed to consider only 
	#weighted edges.
	return nodQDgree,nbMiss,nbcMiss,nbConnctin


########################################################
#~ Compute the actual number of neighbor missed nodes. These values are 
# required to intersect networks. 
########################################################
def calcMiss(G_ntwrkxA,G_ntwrkxB,nod):
	"""
	Input: G_ntwrkxA is one graph of interest, for instance generated 
	from the endogenous competition data. G_ntwrkxB is one graph of 
	interest, for instance generated from the expression data, nod is 
	the node of interest to be tested.
	Output: (Weighted) missed neighbors (nbMiss) and neighbors 
	connections (nbcMiss) in G_ntwrkxB for nod.
	"""
	nbMiss=0.0
	lNeighbornQ=G_ntwrkxA[nod].keys()
	for nodCrrltd in set(lNeighbornQ).difference(set(G_ntwrkxB[nod]. \
	keys())):
		wght=G_ntwrkxA[nod][nodCrrltd]['weight']
		nbMiss+=wght
	#
	nbcMiss=0.0
	for nodNeigbrQA,nodNeigbrQB in combinations(lNeighbornQ,2):
		try:
			wght=G_ntwrkxA[nodNeigbrQA][nodNeigbrQB]['weight']
			try:
				wght2=G_ntwrkxB[nodNeigbrQA][nodNeigbrQB]['weight']
			except:
				nbcMiss+=wght
		except:
			pass
	#
	return nbMiss,nbcMiss


########################################################
#~ Return a set of nodes in two large graphs that share aproximately 
# the same structure (neighbors of neighbors).
########################################################
def growMatch(G_ntwrkxA,G_ntwrkxB,Mimp):
	"""
	Input: G_ntwrkxA is one graph of interest, for instance generated 
	from the endogenous competition data, G_ntwrkxB is one graph of 
	interest, for instance generated from the expression data, Mimp 
	contains the matches for the important nodes in G_ntwrkxA and 
	G_ntwrkxB.
	Output: M contains the node matches for the resulting graph match
	put all node matches from Mimp to a priority queue Q sorted by their
	qualities.
	"""
	#----------------------------
	#~ Core Method to match nodes
	def matchNodes(Sq,Sdb,Q):
		"""
		Input:Sq is a set of nodes in G_ntwrkxA, Sdb is a set of nodes 
		in G_ntwrkxB, M contains all the current node matches found so 
		far, Q contains all the candidate node matches to be examined.
		"""
		for nodQ in Sq:
			if nodQ in Sdb:
				if nodQ not in Q:
					Q.append(nodQ)
				Sdb.remove(nodQ)
	#----------------------------
	#~ Calculate nodes that share neighbors and neighbors of neighbors
	def examineNodesNearBy(G_ntwrkxA,G_ntwrkxB,bstQualNode,M,Q, 
		exact=True):
		"""
		Input: G_ntwrkxA is the graph generated from the endogenous 
		competition data, G_ntwrkxB is the graph generated from the 
		expression data, bstQualNode is a node in G_ntwrkxA and 
		G_ntwrkxB, M contains all the current node matches found so 
		far, Q contains all the candidate node matches to be examined.
		Note: If the variable exact is off nodes connected up to one 
		neighborn are returned.
		"""		
		NB1q=set(G_ntwrkxA[bstQualNode].keys()).difference(M)#immediate 
		#neighbors of bstQualNode in G_ntwrkxA thathaveno matches in Mc
		NB1db=set(G_ntwrkxB[bstQualNode].keys()).difference(M.union(Q))
		#immediate neighbors of bstQualNode in G_ntwrkxB that have no 
		#matches in Mc and Q
		matchNodes(NB1q,NB1db,Q)
		if not exact:#look for neighbors of neighbors
			NB2q=set()#nodes two hops away from Nq that haveno matches 
			#in Mc
			for neighbNod in NB1q:
				NB2q.update(set(G_ntwrkxA[neighbNod].keys()))
			NB2q=NB2q.difference(M)
			#
			NB2db=set()#nodes two hopsaway from Nq that haveno matches 
			#in Mc
			for neighbNod in NB1db:
				NB2db.update(set(G_ntwrkxB[neighbNod].keys()))
			NB2db=NB2db.difference(M.union(Q))
			#
			matchNodes(NB1q,NB2db,Q)# Silenced in case of exact matches
			matchNodes(NB2q,NB1db,Q)# Silenced in case of exact matches
	#----------------------------
	#put all nodes in Mimp to a priority queue Q sorted by their 
	#qualities
	Q=sorted([(nod.w,nod.name) for nod in Mimp])#only includes a small 
	#set of important nodes
	Q=[nod[1] for nod in Q]
	M=set()
	while Q:
		bstQualNode=Q.pop()
		M.add(bstQualNode)
		examineNodesNearBy(G_ntwrkxA,G_ntwrkxB,bstQualNode,M,Q)
	return M


########################################################
#~ Return two subgraphs sharing a similar structure (neigbors of 
# neigbors) from two networks of interest.
########################################################
def rtrnSbgrphsApprxStrctr(G_ntwrkxA,G_ntwrkxB,rho,theta,nThrd, \
	vrbse=True):
	"""
	Input: G_ntwrkxA is one graph of interest, for instance generated 
	from the endogenous competition data, G_ntwrkxB is one graph of 
	interest, for instance generated from the expression data, rho is 
	the degree of approximation orthe percentage of neighbors of a query 
	node that can have no corresponding matches in the neighborhood of a 
	database node, and theta is equivalent to rho, just that it is 
	applicable to the sum of the weighted edges between neighbors. 
	Neighbor connections are the number of edges between the neighbors. 
	If vrbse is True, the log forthe run is going to be printed.
	Output: G_ntwrkxSbgrphA is the subgraph of G_ntwrkxA with important 
	nodesand weighted edges connecting them in G_ntwrkxB,G_ntwrkxSbgrphB 
	is the subgraph of G_ntwrkxB with the important nodes and weighted 
	edges in G_ntwrkxA. ovrAllogRun is the log of the run.
	Note: The edges inthe output are only those shared by the two graphs
	between the important nodes.
	"""
	#----------------------------
	# return all edges shared by important nodes
	def rtrnAllEdgs(G_ntwrkxA,G_ntwrkxB,M):
		"""
		Input: G_ntwrkxA is the graph generated from the endogenous 
		competition data, G_ntwrkxB is the graph generated from the 
		expression data, Sq is a set of nodes in G_ntwrkxA, Sdb is a 
		set of nodes in G_ntwrkxB, M contains all the current node 
		matches found.
		Output: G_ntwrkxA_ovrlp is the subgraph of G_ntwrkxA with the
		nodes in M, G_ntwrkxB_ovrlp is the subgraph of G_ntwrkxB 
		with the nodes in M.
		Note: Only the edges shared by the nodes in M in the two graphs 
		are returned.
		"""
		sG_ntwrkxA_nds=set(G_ntwrkxA.nodes())
		sG_ntwrkxB_nds=set(G_ntwrkxB.nodes())
		del_sG_ntwrkxA_nds=sG_ntwrkxA_nds.difference(set(M))
		del_sG_ntwrkxB_nds=sG_ntwrkxB_nds.difference(set(M))
		G_ntwrkxA.remove_nodes_from(del_sG_ntwrkxA_nds)
		G_ntwrkxB.remove_nodes_from(del_sG_ntwrkxB_nds)
		for nodA,nodB in combinations(M,2):
			if G_ntwrkxA.has_edge(nodA,nodB) and not \
			G_ntwrkxB.has_edge(nodA,nodB):
				G_ntwrkxA.remove_edge(nodA,nodB)
			elif G_ntwrkxB.has_edge(nodA,nodB) and not \
			G_ntwrkxA.has_edge(nodA,nodB):
				G_ntwrkxB.remove_edge(nodA,nodB)
		return G_ntwrkxA,G_ntwrkxB
	#----------------------------
	# return nodes in subgraphs sharing a similar structure
	Mimp,ovrAllogRun = slctImprtntNds(G_ntwrkxA,G_ntwrkxB,rho,theta, \
	nThrd,vrbse)
	M = growMatch(G_ntwrkxA,G_ntwrkxB,Mimp)
	#----------------------------
	# report results on M
	try:
		raise exceptions.CelleryWarningObjct \
		('%s important nodes (including neighbors)'%len(M), \
		'were found to be shared by the two graphss')
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	#----------------------------
	# return subgraphs sharing a similar structure (i.e. nodes in M)
	G_ntwrkxSbgrphA,G_ntwrkxSbgrphB = rtrnAllEdgs(G_ntwrkxA,G_ntwrkxB,M)
	return G_ntwrkxSbgrphA,G_ntwrkxSbgrphB,ovrAllogRun


########################################################
#~ Find a minimum (maximum given weights) spaning tree for a list of
# edges
########################################################
def mst_krukal(sNds,lWghtEdgs,vrbse=True):
	"""
	Input: sNds is a set of nodes in a graph of interest, lWghtEdgs is
	a list of of tuples in the form [(w0,edg0),(w1,edg1),...(wN,edgN)],
	where wX is the weight of edgeX and edgeX has the structure ndA,ndB
	and {ndA,ndB} in sNds.If vrbse is True, the log for the run is going 
	to be printed.
	Output: G is a graph in networkX format with the minimum spaning 
	tree built from the input nodes and edges. ovrAllogRun is the log of 
	the run.
	NOTE: Given weights, the algorithm proceeds to maximize the edge
	weigth instead of minimize it.
	"""
	#
	G=nx.Graph()#~ A = []#holder for the valid edges in the MST
	lDisjntSet = __class__.disjntSet(sNds)#make set for each node
	lWghtEdgs.sort()#sort edges in a decreasing order by weight
	lWghtEdgs.reverse()
	lenlWghtEdgs = len(lWghtEdgs)
	cnt = -1
	for wX,ndA,ndB in lWghtEdgs:
		cnt+=1
		if cnt%10000==0:
			print cnt,lenlWghtEdgs
		if lDisjntSet.findSet(ndA)!=lDisjntSet.findSet(ndB):
			G.add_edge(ndA,ndB,{'weight':str(wX)})#A.append((ndA,ndB,wX))
			lDisjntSet.djntUnion(ndA,ndB)
	#----------------------------
	# report results
	ovrAllogRun = []
	try:
		raise exceptions.CelleryWarningObjct \
		('A minimum spanning tree was build with %s nodes',' and %s edges())
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return G,ovrAllogRun


########################################################
#~ Retrieve triangles between edges whose relation can be explain by 
# a third value
########################################################
def rtrnTrngls(G_ntwrkxA,G_ntwrkxB,G_ntwrkxC):
	"""
	Input: G_ntwrkxA is an undirected graph (for instance of ECP 
	-Endogenous Competition Potential-) with important nodes and 
	weighted edges connecting them in G_ntwrkxB (and viceversa). 
	G_ntwrkxC includes all the edges with expression values, 
	independently if they are shared with G_ntwrkxA or G_ntwrkxB.
	Output: G_ntwrkxATrngls is G_ntwrkxA graph with edges between 
	lncRNA-lncRNA and gene-gene present in G_ntwrkxC, this as long 
	as there is lncRNA-gene relation in G_ntwrkxA, G_ntwrkxBTrngls is
	the same for G_ntwrkxB.
	"""
	G_ntwrkxATrngls=G_ntwrkxA.copy()
	G_ntwrkxBTrngls=G_ntwrkxB.copy()
	for nd in G_ntwrkxA.nodes():#same nodes as in G_ntwrkxB
		lNeigbrNds=G_ntwrkxA[nd].keys()
		for nodNeigbrA,nodNeigbrB in combinations(lNeigbrNds,2):
			if G_ntwrkxC.has_edge(nodNeigbrA,nodNeigbrB) and not \
			G_ntwrkxA.has_edge(nodNeigbrA,nodNeigbrB):
				G_ntwrkxATrngls.add_edge(nodNeigbrA,nodNeigbrB, \
				{'weight':G_ntwrkxC[nodNeigbrA][nodNeigbrB] \
				['weight'],'color':'red'})
				G_ntwrkxBTrngls.add_edge(nodNeigbrA,nodNeigbrB, \
				{'weight':G_ntwrkxC[nodNeigbrA][nodNeigbrB] \
				['weight'],'color':'red'})
	#----------------------------
	# complete color edges
	for nodA,nodB in G_ntwrkxATrngls.edges():
		try:
			color=G_ntwrkxATrngls[nodA][nodB]['color']
		except:
			G_ntwrkxATrngls[nodA][nodB]['color']='black'
			G_ntwrkxBTrngls[nodA][nodB]['color']='black'
	return G_ntwrkxATrngls,G_ntwrkxBTrngls


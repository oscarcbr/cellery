#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  stats.py part of cellery (ceRNAs linking inference)
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
Libraries to run statistical testing
"""

########################################################
#~ Import libraries.
########################################################
from networkx.algorithms import bipartite
from numpy import array,float32,isnan,nan,zeros


#----------------------------
#~ p-value correction methods
#----------------------------

########################################################
#~ Correct p-values using Benjamini-Hochberg FDR approach 
########################################################
def bh_qvalues(pVls):
	"""
	Input: pVls is a list or array of p-values with no correction.
	Output: qVls is an array or list of Benjamini-Hochberg FDR q-values 
	corresponding to p-values.
	NOTE: This method allows for nan p-values as input, ignores them
	for calculations, and returns q-values of nan for them.
	NOTE: Output results follow the same datatype as the input.
	"""
	#----------------------------
	#Set initial variables
	if not isinstance(pVls,list):#is an array
		pVls = pVls.tolist()
	qVls = zeros(len(pVls),dtype=float32)
	qVls.fill(nan)
	pValsNonNan,absltPos = zip(*sorted([(val,pos) for pos,val in \
	enumerate(pVls) if not isnan(val)]))
	if pValsNonNan[0] < 0 or pValsNonNan[-1] > 1:
		raise ValueError("p-values must be between 0 and 1")
	m = len(pValsNonNan)
	#----------------------------
	#Run calculations
	mincoeff = pValsNonNan[-1]
	qVls[absltPos[-1]] = mincoeff
	for j in xrange(m-2, -1, -1):
		coeff = (m*pValsNonNan[j])/(j+1)
		if coeff < mincoeff:
			mincoeff = coeff
		qVls[absltPos[j]] = mincoeff
	if isinstance(pVls,list):#is an array
		return qVls.tolist()
	else:
		return qVls


#----------------------------
#~ Method to estimate significance of network structures
#----------------------------	

########################################################
#~ Naive speed algorithm to randomize networks.
########################################################
def naiveSpeed(G_target,n_wrmngStps,n_rndmStps,sNdsIntrst=False, \
	smpl_genrt=1000):
	"""
	Input: G_target is the graph of interest, n_rndmStps is the number
	of walk steps to run randomizations, n_wrmngStps is the number of
	warming steps after which random graphs n_rndmStps-n_wrmngStps are
	going to be recorded.
	Output: lRndmzdG is a list of randomized graphs with the same degree 
	sequences as G_target.
	NOTE: sNdsIntrst is a dummy parameter and always is False.
	"""
	lRndmzdG = []
	G_rndmzd = G_target.copy()
	wrmngTrhsld = n_rndmStps-n_wrmngStps
	gnCnt = smpl_genrt-1#to sample as soon as reach the warming point
	while n_rndmStps>0:
		edges = G_rndmzd.edges()
		edgeIJ,edgeKL = choice(edges),choice(edges)
		if edgeIJ!=edgeKL:
			nodI,nodJ = edgeIJ
			nodK,nodL = edgeKL
			if not G_rndmzd.has_edge(nodI,nodL) and \
			not G_rndmzd.has_edge(nodK,nodJ):
				G_rndmzd.add_edge(nodI,nodL)
				G_rndmzd.add_edge(nodK,nodJ)
				G_rndmzd.remove_edge(nodI,nodJ)
				G_rndmzd.remove_edge(nodK,nodL)
				#
				n_rndmStps-=1
				if n_rndmStps%1000==0:
					print 'running warming step %s'%n_rndmStps
				#Save graph after convergence == warming steps
				if n_rndmStps<wrmngTrhsld:
					gnCnt+=1
					if gnCnt==smpl_genrt:
						lRndmzdG.append(G_rndmzd)
						gnCnt = 0
	return lRndmzdG


########################################################
#~ self_loop algorithm to randomize networks.
########################################################
def self_loop(G_target,n_wrmngStps,n_rndmStps,sNdsIntrst=False, \
	smpl_genrt=1000):
	"""
	Input: G_target is the graph of interest, n_rndmStps is the number
	of walk steps to run randomizations, n_wrmngStps is the number of
	warming steps after which random graphs n_rndmStps-n_wrmngStps are
	going to be recorded. The value smpl_genr determines after the 
	warming steps how often the samples are to be taken.
	Output: lRndmzdG is a list of randomized graphs with the same degree 
	sequences as G_target following a Markov-chain approach as non-
	swapable pairs are taken as self-loops. This are expected to be 
	optimal after n_wrmngStps steps.
	Note: this method was modified to fit the self_loop algorithm in 
	Gionis et al. by using the method double_edge_swap in networkX.
	Needed to add a smpl_genrt value as all the samples were the same.
	NOTE: sNdsIntrst is a dummy parameter and always is False.
	"""
	lRndmzdG = []
	G_rndmzd = G_target.copy()
	wrmngTrhsld = n_rndmStps-n_wrmngStps
	gnCnt = smpl_genrt-1#to sample as soon as reach the warming point
	while n_rndmStps>0:
		n_rndmStps-=1
		if n_rndmStps%1000==0:
			print 'running randomization step %s'%n_rndmStps
		edges = G_rndmzd.edges()
		edgeIJ,edgeKL = choice(edges),choice(edges)
		if edgeIJ!=edgeKL:
			nodI,nodJ = edgeIJ
			nodK,nodL = edgeKL
			if not G_rndmzd.has_edge(nodI,nodL) and \
			not G_rndmzd.has_edge(nodK,nodJ):
				G_rndmzd.add_edge(nodI,nodL)
				G_rndmzd.add_edge(nodK,nodJ)
				G_rndmzd.remove_edge(nodI,nodJ)
				G_rndmzd.remove_edge(nodK,nodL)
		#Save graph after convergence == warming steps
		if n_rndmStps<wrmngTrhsld:
			gnCnt+=1
			if gnCnt==smpl_genrt:
				lRndmzdG.append(G_rndmzd)
				gnCnt = 0
	return lRndmzdG


########################################################
#~ metropolis_hastings algorithm to randomize networks.
########################################################
def metropolis_hastings(G_target,n_wrmngStps,n_rndmStps,sNdsIntrst, \
	smpl_genrt=100):
	"""
	Input: G_target is the graph of interest, n_rndmStps is the number
	of walk steps to run randomizations, n_wrmngStps is the number of
	warming steps after which random graphs n_rndmStps-n_wrmngStps are
	going to be recorded, sNdsIntrst is the set of "left" nodes in the 
	bipartite graph.
	Output: lRndmzdG is a list of randomized graphs with the same degree 
	sequences as G_target following a Markov-chain with a Metropolis-
	Hastings approach. This are expected to be optimal after n_wrmngStps 
	steps.
	NOTE: This method requires a bipartite graph. To ensure the graph is
	bipartite, a graph requires a set of nodes of interest. This willfor
	sure work in case were there is lncRNA-gene connections, but should
	be tested case by case if lncRNA-lncRNA gene-gene connections are 
	included.
	"""
	assert bipartite.is_bipartite(G_target)
	#----------------------------
	#~ Find adjacent method
	def find_adjacent(G_rndmzd):
		"""
		Input: G_rndmzd is the graph of interest in process of 
		randomization
		Output: G_rndmzd_prime is a graph that differs from G_rndmzd in 
		exactly one swap (i.e., (G,G') in T).
		"""
		G_rndmzd_prime = G_rndmzd.copy()
		edges = G_rndmzd.edges()
		while True:
			edgeIJ,edgeKL = choice(edges),choice(edges)
			if edgeIJ! = edgeKL:#works well as edges are sorte already in 
				#"edges"
				nodI,nodJ = edgeIJ
				nodK,nodL = edgeKL
				if not G_rndmzd.has_edge(nodI,nodL) and \
				not G_rndmzd.has_edge(nodK,nodJ):
					G_rndmzd_prime.add_edge(nodI,nodL)
					G_rndmzd_prime.add_edge(nodK,nodJ)
					G_rndmzd_prime.remove_edge(nodI,nodJ)
					G_rndmzd_prime.remove_edge(nodK,nodL)
					break	
		return G_rndmzd_prime
	#----------------------------
	#~ Run	
	lRndmzdG = []
	G_rndmzd = G_target.copy()
	dG_rndmzd = abs(clcDgrph(G_rndmzd,sNdsIntrst))
	wrmngTrhsld = n_rndmStps-n_wrmngStps
	gnCnt = smpl_genrt-1#to sample as soon as reach the warming point
	while n_rndmStps>0:
		n_rndmStps-=1
		G_rndmzd_prime = find_adjacent(G_rndmzd)
		dG_rndmzd_prime = abs(clcDgrph(G_rndmzd_prime,sNdsIntrst))
		acceptncProb = min(1,(dG_rndmzd/dG_rndmzd_prime))
		if n_rndmStps%1000==0:
			print 'running randomization step %s'%n_rndmStps
		if bernoulli(acceptncProb).rvs():#returns 1 (acceptance) 
			#following a bernoulli distribution with p=acceptncProb
			G_rndmzd = G_rndmzd_prime
			dG_rndmzd = abs(clcDgrph(G_rndmzd,sNdsIntrst))
		#Save graph after convergence == warming steps
		if n_rndmStps<wrmngTrhsld:
			gnCnt+=1
			if gnCnt==smpl_genrt:
				lRndmzdG.append(G_rndmzd)
				gnCnt=0
	return lRndmzdG


########################################################
#~ clcDgrph method for metropolis_hastings algorithm to randomize 
# networks.
########################################################
def clcDgrph(G,sNdsIntrst):
	"""
	Input: G is a graph of interest MUST be bipartite.
	Output: dG is the number of swappable pairs of graph G.
	Note: dG is calculated as Ginois et al. 2007.
	"""
	#----------------------------
	#Given a bipartite graph G = (U,V,E)
	def J(G,row_order,column_order):
		"""
		Input: G is a graph of interest, row_order is the list of left
		nodes and column_order is the list of right nodes in a bipartite
		graph.
		Output: J_G is the number of disjoint pairs of edges.
		"""
		lenEdgs = G.number_of_edges()
		sumR2 = sum([G.degree(i)**2 for i in row_order])
		sumC2 = sum([G.degree(j)**2 for j in column_order])
		J_G = 0.5*((lenEdgs*(lenEdgs+1)) - sumR2 - sumC2)
		return J_G
	#----------------------------
	def Z(G,row_order,column_order):
		"""
		Input: G is a graph of interest, row_order is the list of left
		nodes and column_order is the list of right nodes in a bipartite
		graph.
		Output: Z_G is the number of "Z" structures.
		"""
		Z_G=0.0
		for i in row_order:
			riM1=G.degree(i)-1
			for j in G[i].keys():#assume row are only connected to 
				#columns
				cjM1=G.degree(j)-1
				Z_G+=riM1*cjM1
		return Z_G
	#----------------------------
	def K22(M,row_order,row):
		"""
		Input: M = DD^T is a graph of interest, row_order is the list of
		left nodes.
		Output: K22_G is the number of cliques of G
		"""
		K22_G=0.0
		for i in row_order:
			for k in row_order:
				if i!=k:
					M_ik=M[row[i],row[k]]
					K22_G+=(M_ik**2)-M_ik
		#K22_G*=0.5 silenced to avoid the multiplication by 2 when 
		#calculating d(G)
		return K22_G
	#----------------------------
	#Calculate dG now
	row_order=sorted(sNdsIntrst)
	column_order = list(set(G) - sNdsIntrst)
	row = dict(zip(row_order,count()))
	col = dict(zip(column_order,count()))
	D = np.zeros((len(row),len(col)))
	for u in row_order:
		for v in column_order:
			if G.has_edge(u,v):
				D[row[u],col[v]] = 1
	#
	M=D.dot(D.transpose())#DD^T
	J_G=J(G,row_order,column_order)
	Z_G=Z(G,row_order,column_order)
	K22_G=K22(M,row_order,row)
	dG=J_G-Z_G+K22_G
	return dG	


########################################################
#~ Estimate significance of network structures (percentage in which the
# nodes in the target graph (G) appears in the list of randomized
# graphs)
########################################################
def clcStsts(G,lRndmzdG):
	"""
	Input:  G is the graph of interest, lRndmzdG is a list of randomized 
	graphs.
	Output: dEdgesProb is the proportion of times that each edge in G 
	is found in the graphs of lRndmzdG, pGrph is the proportion of 
	graphs in lRndmzdG that are identical to G.
	"""
	pGrph=0
	dEdgesProb=dict([(edge,0) for edge in G.edges()])
	for G_rmdzd in lRndmzdG:
		sameTst=True
		for edge in G.edges():
			if G_rmdzd.has_edge(*edge):
				dEdgesProb[edge]+=1
			else:
				sameTst=False
		if sameTst:
			pGrph+=1
	#
	lenlRndmzdG=float32(len(lRndmzdG))
	dEdgesProb=dict([(edge,float32(dEdgesProb[edge]/lenlRndmzdG)) for \
	edge in dEdgesProb.keys()])#turn counts into frequencies
	pGrph/=lenlRndmzdG
	#
	return dEdgesProb,pGrph


########################################################
#~ Wrap the randomization of a graph and estimate the significance of 
# its structure
########################################################
def wrpRndmzNtwrk(G,algrthm='naiveSpeed',n_wrmngStps=1000, \
	n_rndmStps=100000,sNdsIntrst=False,smpl_genrt=100,vrbse=True):
	"""
	Input: G is the graph of interest. algrthm is the algorithm to 
	randomize the input graph: {'naiveSpeed','self_loop',
	'metropolis_hastings'}. n_rndmStps is the number of walk steps to 
	run randomizations. n_wrmngStps is the number of warming steps after 
	which random graphs n_rndmStps-n_wrmngStps are going to be recorded, 
	sNdsIntrst is the set of "left" nodes in the bipartite graph (only 
	for 'metropolis_hastings'). If vrbse is True, the log for the run is 
	going to be printed.
	Output: ovrAllogRun is the log of the run.
	NOTE: The method runs a randomization over G following the input 
	parameters and annotates G with the resulting probabilities under
	the value for algrthm. Also prints the probability of the network
	structure.
	"""
	ovrAllogRun = []
	lRndmzdG = getattr(statistics,algrthm)(G,algrthm,n_wrmngStps, \
	n_rndmStps,sNdsIntrst,smpl_genrt)
	dEdgesProb,pGrph = statistics.clcStsts(G,lRndmzdG)
	ovrAllogRun_annt = anntNtwrkXEdgs(G_target,dEdgesProb,algrthm)
	ovrAllogRun.extend(ovrAllogRun_annt)
	#----------------------------
	# report final results
	try:
		raise exceptions.CelleryWarningObjct \
		('The probability for the input graph to happen by chance is %s'% \
		pGrph,'as calculated by the algorithm: %s, with parameters: %s.'% \
		(algrthm,'.'.join(['n_wrmngStps: ',n_wrmngStps,'n_rndmStps: ', \
		n_rndmStps,'smpl_genrt: ',smpl_genrt])))
	except exceptions.CelleryWarningObjct as mssge:
		if vrbse:
			print mssge
		ovrAllogRun.append(mssge)
	return ovrAllogRun
	

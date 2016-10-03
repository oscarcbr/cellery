#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_subgraphs.py part of cellery (ceRNAs linking inference)
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


########################################################
#~ Import libraries.
########################################################
from cellery import __object__
from cellery import alignments
import os

########################################################
#~ Construct objects for sequences of protein coding genes, 3'UTRs and 
# lncRNAs from sql databases.
########################################################
def anntGrph(infGns,infTKeg,G_targetGML,anntGO=True):
	"""
	Input: infGns is the extended annotation file from a UCSC table for
	ENSEMBL genes and transcripts with added GO terms, and KEGG pathways
	fields. infTKeg is a KEGG downloaded file with the ENSEMBL gene and
	the KEGG pathways (must be obtained from Galaxy Genome Diversity's
	respositoire). G_targetGML is the graph of interest to annotate in 
	'.gml' format.
	Output: G_annttdGML is going to be written in the file infGns.
	Note: Optionally anntGO can be turned off as it takes a lot of 
	resource to annotate genes with it.
	"""

	#
	#Give me the annotation
	#
	dEnsmbTKEGG=mkDctAnnt(infTKeg,0,2)
	dEnsmbNam=mkDctAnnt(infGns,12,16)
	if anntGO:
		dEnsmbGO=mkDctAnnt(infGns,12,20)
		sGO=prsAnnt(dEnsmbGO)
	#
	sEnsmbl=set(dEnsmbTKEGG.keys())
	#dict of Ensembl gene names to KEGG pathways
	dEnsmbKEGG=dict([(x.split('\t')[12],dEnsmbTKEGG[x.split('\t')[1]]) \
	for x in open(infGns).read().splitlines() if x.strip() and \
	x.split('\t')[1] in sEnsmbl and dEnsmbTKEGG[x.split('\t')[1]]!='N'])
	

	#~ 
	parser = argparse.ArgumentParser(description='This code search for shared edges between two weighted undirected graphs, or alternatively a suboptimal subgraph shared between two graphs. In addition, it provides randomization and biological functional annotation. It takes as input several files with matrices representing expression and competition, but it can be further extended to test other propierties.')
	parser.add_argument('--outf',metavar='output folder',type=str,help='the path to output folder.',required=True)
	parser.add_argument('--ensemblAnntn',metavar='input txt file',type=str,help='the extended annotation file. Can be obtained from the Galaxy repository.',required=False,default=False)
	parser.add_argument('--KEGGAnntn',metavar='input txt file',type=str,help='the KEGG annotation file, it has in one column the KEGG gene code and in the other the KEGG pathway. Can be obtained from the Galaxy repository or from the KEGG web page',required=False,default=False)
	parser.add_argument('--fECPMtrxTrgt',metavar='input ecp file',type=str,help='the endogenous competition potential matrix, it can be generated with the mirna_trgtvsbckGrps code.',required=True)
	parser.add_argument('--fExprCrrRMtrxTrgt',metavar='input crr file',type=str,help='the expression correlation matrix with R values in the target dataset, it can be generated with the mirna_trgtvsbckGrps code.',required=True)
	parser.add_argument('--fExprCrrPMtrxTrgt',metavar='input crp file',type=str,help='the expression correlation matrix with p/q values in the target dataset, it can be generated with the mirna_trgtvsbckGrps code.',required=True)
	parser.add_argument('--fECPMtrxRef',metavar='input ecp file',type=str,help='the endogenous competition potential matrix in the reference dataset, it can be generated with the mirna_trgtvsbckGrps code.',required=False,default=False)
	parser.add_argument('--fExprCrrRMtrxRef',metavar='input crr file',type=str,help='the expression correlation matrix with R values in the reference dataset, it can be generated with the mirna_trgtvsbckGrps code.',required=False,default=False)
	parser.add_argument('--fExprCrrPMtrxRef',metavar='input crp file',type=str,help='the expression correlation matrix with p/q values in the reference dataset, it can be generated with the mirna_trgtvsbckGrps code.',required=False,default=False)
	parser.add_argument('--fldrBckGrnd',metavar='input path to folder with background files',type=str,help='input path to folder with background ecp, crr and crp files.',required=False,default=False)
	parser.add_argument('--rndmz_graph',metavar='boolean randomize target graph',type=str,help='boolean randomize target graph. By default True.',required=False,default='True')
	parser.add_argument('--self_loop',metavar='boolean to use self_loop algorithm',type=str,help='boolean to use self_loop algorithm. By default False.',required=False,default='False')
	parser.add_argument('--naive',metavar='boolean to use naive algorithm',type=str,help='boolean to use naive algorithm. By default False',required=False,default='False')
	parser.add_argument('--metropolis_hastings',metavar='boolean to use metropolis_hastings algorithm',type=str,help='boolean to use naive algorithm. By default True.',required=False,default='True')
	parser.add_argument('--rho',metavar='float tolerance to missed neighbor nodes',type=float,help='tolerance to missed neighbor nodes. 0-1 where 1 is more tolerance.',required=False,default=False)
	parser.add_argument('--theta',metavar='float tolerance to missed neighbor weighted edges',type=float,help='tolerance to missed neighbor weighted edges. 0-1 where 1 is more tolerance.',required=False,default=False)
	parser.add_argument('--sgnTrhld',metavar='float significance to trim probabilities in expressions',type=float,help='significance to trim probabilities in expressions.',required=True)
	parser.add_argument('--prctECPcut',metavar='int bottom percentile to cut Endogenous Competitions Potential',type=int,help='bottom percentile to cut Endogenous Competitions Potential values (0-100).',required=True)
	parser.add_argument('--prctExprCrrCut',metavar='int bottom percentile to cut expression values',type=int,help='bottom percentile to cut expression values (0-100).',required=True)
	parser.add_argument('--n_rndmStps_aftrWrmng',metavar='int number of graphs to be safe after heating steps',type=int,help='number of graphs to be analyze after heating steps. By default 1000.',required=False,default=100000)
	parser.add_argument('--smpl_genrt',metavar='int number of generations to take one sample',type=int,help='interval in generations to sample. By default 100.',required=False,default=100)
	parser.add_argument('--n_wrmngStps',metavar='int number of warming steps to allow graph to mix',type=int,help='number of warming steps to allow graph to mix. By default determine as 6 times the number of edges in target graph.',required=False,default=False)
	parser.add_argument('--nThrd',metavar='int number of threads to run the subgraph algorithm',type=int,help='int number of threads to run the subgraph algorithm. Default 5',required=False,default=5)
	parser.add_argument('--optmlSbgrph',metavar='bool use the suboptimal subgraph algorithm from Tian 2008',type=str,help='bool use the suboptimal subgraph algorithm from Tian 2008 (i.e. TALE). If false, the subgraph is going to be selecting by choosing only those edges present in both graphs and then deleting all the nodes wiht 0-degree. By default',required=False,default='False')
	parser.add_argument('--fGMLTrnglsECP',metavar='output gml file with triangles in ECP (not tested)',type=str,help='this includes the subgraph with edges of triangles colored in red when supported by the expression graph. It is optional and should consider all vs. all matrices. By default False.',required=False,default=False)
	parser.add_argument('--fGMLTrnglsExprCrr',metavar='output gml file with triangles in expression values (not tested)',type=str,help='this includes the subgraph with edges of triangles colored in red when supported by the ECP graph. It is optional and should consider all vs. all matrices. By default False.',required=False,default=False)
	args = parser.parse_args()
	#~ 
	outfldr = args.outf#Output folder
	#Annotation files
	infGns = args.ensemblAnntn#annotation table
	infTKeg = args.KEGGAnntn#KEGG annotation
	#only files having the three extension '.crr', '.crp' and '.ecp' are 
	#going to be taken into account.
	#Target files
	fECPMtrxTrgt = args.fECPMtrxTrgt
	fExprCrrRMtrxTrgt = args.fExprCrrRMtrxTrgt
	fExprCrrPMtrxTrgt = args.fExprCrrPMtrxTrgt
	#Reference files
	fECPMtrxRef = args.fECPMtrxRef
	fExprCrrRMtrxRef = args.fExprCrrRMtrxRef
	fExprCrrPMtrxRef = args.fExprCrrPMtrxRef
	fldrBckGrnd = args.fldrBckGrnd
	assert fldrBckGrnd or (fECPMtrxRef and fExprCrrRMtrxRef and fExprCrrPMtrxRef)
	#Randomize target graph
	rndmz_graph = args.rndmz_graph
	if upper(rndmz_graph)=='FALSE':
		rndmz_graph=False
	elif upper(rndmz_graph)=='TRUE':
		rndmz_graph=True
	alSelf_loop = args.self_loop
	if upper(alSelf_loop)=='FALSE':
		alSelf_loop=False
	elif upper(alSelf_loop)=='TRUE':
		alSelf_loop=True
	alNaive = args.naive
	if upper(alNaive)=='FALSE':
		alNaive=False
	elif upper(alNaive)=='TRUE':
		alNaive=True
	alMetropolis_hastings = args.metropolis_hastings
	if upper(alMetropolis_hastings)=='FALSE':
		alMetropolis_hastings=False
	elif upper(alMetropolis_hastings)=='TRUE':
		alMetropolis_hastings=True
	if rndmz_graph:
		assert alSelf_loop or alNaive or alMetropolis_hastings
	#Parameters
	rho = args.rho
	theta = args.theta
	sgnTrhld = args.sgnTrhld
	prctECPcut = args.prctECPcut
	prctExprCrrCut = args.prctExprCrrCut
	n_rndmStps_aftrWrmng = args.n_rndmStps_aftrWrmng
	smpl_genrt = args.smpl_genrt
	n_wrmngStps = args.n_wrmngStps
	assert 0<=rho<=1 and 0<=theta<=1 and 0<=sgnTrhld<=1
	assert 0<=prctExprCrrCut<=100 and 0<=prctECPcut<=100
	#Output graph files
	f = os.path.split(fECPMtrxTrgt)[1].split('.')[0]
	fExprCrrGrphTrgt = os.path.join(outfldr,'%s_crr.gml'%f)
	fECPGrphTrgt = os.path.join(outfldr,'%s_ecp.gml'%f)
	#this includes the subgraph with edges of triangles colored in red 
	#when supported by the expression graph. It is optional and should 
	#consider all vs. all matrices
	fGMLTrnglsECP = args.fGMLTrnglsECP
	fGMLTrnglsExprCrr = args.fGMLTrnglsExprCrr
	nThrd=args.nThrd
	#Select algorithm to substract subgraph
	optmlSbgrph = args.optmlSbgrph
	if upper(optmlSbgrph)=='FALSE':
		optmlSbgrph=False
	elif upper(optmlSbgrph)=='TRUE':
		optmlSbgrph=True
	if optmlSbgrph:
		assert rho and theta
	#Target - suboptimal subgraph between ECP and expression
	Gendo_compTrgt,Gcorrl_exprTrgt,sLncRNAallTrgt = wrprMkGrph(fECPMtrxTrgt, \
	fExprCrrRMtrxTrgt,fExprCrrPMtrxTrgt,fECPGrphTrgt,fExprCrrGrphTrgt, \
	prctECPcut,prctExprCrrCut,rho,theta,sgnTrhld,optmlSbgrph,nThrd, \
	fGMLTrnglsECP,fGMLTrnglsExprCrr)
	#~ 
	if not n_wrmngStps:
		n_wrmngStps = 6*Gendo_compTrgt.number_of_edges()#5L swaps is the step 
		#in which convergence initiates at latest -as indicated by Gioinis (L=edges)-
	#~ Process Background
	if fldrBckGrnd:
		lFlsBckGrnd=set([os.path.splitext(f)[0] for f in os.listdir(fldrBckGrnd)])
		lFlsBckGrnd=[f for f in lFlsBckGrnd if os.path.exists(os.path.join \
		(fldrBckGrnd,'%s.crr'%f)) and os.path.exists(os.path.join(fldrBckGrnd, \
		'%s.crp'%f)) and os.path.exists(os.path.join(fldrBckGrnd,'%s.ecp'%f))]
		colorGraph=False#for next coloring step
	else:#The reference files must have the same prefix
		fldrBckGrnd,f=os.path.split(fECPMtrxRef)
		lFlsBckGrnd=[os.path.splitext(f)[0]]
		colorGraph=True
	#Make background
	#In the original case: fExprCrrPMtrxTrgt==fExprCrrPMtrxRef AND 
	#fExprCrrRMtrxTrgt==fExprCrrRMtrxRef
	lBckGrnd_G_ECP=[]#lBckGrnd_G_Crr=[] can be skipped as lBckGrnd_G_ECP and 
	#lBckGrnd_G_Crr has the same edges.
	for f in lFlsBckGrnd:
		fECPMtrxRef = os.path.join(fldrBckGrnd,'%s.ecp'%f)
		fExprCrrRMtrxRef = os.path.join(fldrBckGrnd,'%s.crr'%f)
		fExprCrrPMtrxRef = os.path.join(fldrBckGrnd,'%s.crp'%f)
		fECPGrphRef = os.path.join(outfldr,'%s_ecp.gml'%f)
		fExprCrrGrphRef = os.path.join(outfldr,'%s_crr.gml'%f)
		#Reference - suboptimal subgraph between ECP and expression
		Gendo_compRef,Gcorrl_exprRef,sLncRNAallRef = wrprMkGrph(fECPMtrxRef, \
		fExprCrrRMtrxRef,fExprCrrPMtrxRef,fECPGrphRef,fExprCrrGrphRef, \
		prctECPcut,prctExprCrrCut,rho,theta,sgnTrhld,optmlSbgrph,nThrd)
		#Color differences, only for single background graphs
		if colorGraph:
			#Now color edge differences in reference and target graphs in red
			completeGraph(Gendo_compTrgt,Gendo_compRef)
			completeGraph(Gcorrl_exprTrgt,Gcorrl_exprRef)
		else:
			lBckGrnd_G_ECP.append(Gendo_compRef)
			#lBckGrnd_G_Crr.append(Gcorrl_exprRef) can be skipped
	#
	if lBckGrnd_G_ECP:
		dEdgesProbBck,pGrphBck = clcStsts(Gendo_compTrgt,lBckGrnd_G_ECP)
		anntEdgs(Gendo_compTrgt,dEdgesProbBck,'background')
		anntEdgs(Gcorrl_exprTrgt,dEdgesProbBck,'background')
		
	#At this point is important to mention that the randomizations are
	#going to be conducted on the Gendo_compTrgt graph. Because these
	#are 1-0 edges and not wieghted edges, the graph is the same as in
	#Gcorrl_exprTrgt.
	if rndmz_graph:
		#Randomize graphs
		n_rndmStps = n_wrmngStps+n_rndmStps_aftrWrmng 
		#don't use self_loop, it sucks
		if alSelf_loop:
			lRndmzdG1 = self_loop(Gendo_compTrgt,n_wrmngStps,n_rndmStps,smpl_genrt)
			#Now calculate probabilities of edges and the graph
			dEdgesProb1,pGrph1 = clcStsts(Gendo_compTrgt,lRndmzdG1)
			#Annotate edges with probabilities
			anntEdgs(Gendo_compTrgt,dEdgesProb1,'self_loop')
			anntEdgs(Gcorrl_exprTrgt,dEdgesProb1,'self_loop')
		#don't use naiveSpeed==naive, it sucks as well
		if alNaive:
			lRndmzdG2 = naiveSpeed(Gendo_compTrgt,n_wrmngStps,n_rndmStps,smpl_genrt)
			dEdgesProb2,pGrph2 = clcStsts(Gendo_compTrgt,lRndmzdG2)
			anntEdgs(Gendo_compTrgt,dEdgesProb2,'naive')
			anntEdgs(Gcorrl_exprTrgt,dEdgesProb2,'naive')
		#If the Metropolis-Hasting algorithm is going to be used, it is 
		#required that the graph is bipartited
		if alMetropolis_hastings:
			sLncRNAallTrgt=sLncRNAallTrgt.intersection(set(Gendo_compTrgt.nodes()))
			lRndmzdG3 = metropolis_hastings(Gendo_compTrgt,n_wrmngStps,n_rndmStps, \
			sLncRNAallTrgt,smpl_genrt)
			dEdgesProb3,pGrph3 = clcStsts(Gendo_compTrgt,lRndmzdG3)
			anntEdgs(Gendo_compTrgt,dEdgesProb3,'metropolis_hastings')
			anntEdgs(Gcorrl_exprTrgt,dEdgesProb3,'metropolis_hastings')
	#the names make that the previously built images to be replaced
	nx.write_gml(Gendo_compTrgt,fECPGrphTrgt)
	nx.write_gml(Gcorrl_exprTrgt,fExprCrrGrphTrgt)
	if not lBckGrnd_G_ECP:
		nx.write_gml(Gendo_compRef,fECPGrphRef)
		nx.write_gml(Gcorrl_exprRef,fExprCrrGrphRef)
	#Annotate graphs
	if infGns and infTKeg:
		anntGrph(infGns,infTKeg,fECPGrphTrgt,anntGO=False)
		anntGrph(infGns,infTKeg,fECPGrphRef,anntGO=False)
		anntGrph(infGns,infTKeg,fExprCrrGrphTrgt,anntGO=False)
		anntGrph(infGns,infTKeg,fExprCrrGrphRef,anntGO=False)
	#~ 
	#
	print 'Done!'
	if rndmz_graph:
		print 'Results:'
		if alSelf_loop:
			print '\tGraph probability calculated with "self_loop" algorithm: %s'%pGrph1
		if alNaive:
			print '\tGraph probability calculated with "naive" algorithm: %s'%pGrph2
		if alMetropolis_hastings:
			print '\tGraph probability calculated with "metropolis_hastings" algorithm: %s'%pGrph3
	#
	return 0

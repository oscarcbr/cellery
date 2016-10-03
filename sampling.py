#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sampling.py part of cellery (ceRNAs linking inference)
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
Libraries with model estimation methods
"""


########################################################
#~ Import external libraries
########################################################
#----------------------------
#Import matplotlib libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#----------------------------
#Import other libraries
from cellery import exceptions
from copy import copy
from numpy import absolute,array,divide,flatnonzero,float32,inf, \
isfinite,isnan,linspace,log,ma,min,multiply,nan,warnings,zeros
from scipy import stats
warnings.filterwarnings('ignore')


########################################################
#~ Select model and parameters that maximize the likelihood of data
########################################################
def clcMdl(lLenDt,cntLns,mdl,aprioriParams):
	"""
	Input: lLenDt is a list with the lengths of the data of interest. 
	cntLns are the number of lengths in lLenDt. mdl is the model 
	selected. aprioriParams are the parameters for the aprioriMdl.
	Output: AICc is the Akaike Information Criteria corrected for the 
	log likelihood of the model and parameteres that maximizes the 
	likelihood of data. D is Kolmogorov-Smirnov D value for the model 
	and parameteres that maximizes the likelihood of data. p is the 
	Kolmogorov-Smirnov p value for D. lnLklhd is the log likelihood for 
	tje model and parameteres that maximizes the likelihood of data. 
	params are the parameters maximizes the likelihood of data.
	NOTE: aprioriParams and model shall follow scipy standard names.
	"""
	#----------------------------
	#Fit model and parameters
	if aprioriParams:
		params = aprioriParams
	else:
		params = mdl.fit(lLenDt)
	nparams = len(flatnonzero(array(params)))
	pdfMdl = mdl.pdf(lLenDt,*params)
	lg = log(pdfMdl)
	mask = isfinite(lg)
	lnLklhd = lg[mask].sum()
	#----------------------------
	#Calculate AICc
	AIC = multiply(2,nparams) - multiply(2,lnLklhd)
	if cntLns>nparams+1:
		#akaike_icc
		AICc = divide(AIC+(2*nparams*(nparams+1)),(cntLns-nparams-1))
	else:
		AICc = nan
	D,p = stats.kstest(lLenDt,mdl.name,args=params)
	return AICc,D,p,lnLklhd,params


########################################################
#~ Fit model to the length of a set of objects with a pointer
########################################################	
def fitMdl(lDtObjcts,outPlt=False,pntr=False,sTrgtPntr=False, \
	aprioriMdl=False,aprioriParams=False,vrbose=True):
	"""
	Input:  lDtObjcts is a list with gene/lncrna objects with pointers  
	"len" and the pointer for genes/lncrnas of interest. outPlt is a 
	file name to plot the distribution of lengths and fitted model. pntr 
	is the pointer to the attribute of interest (to select a subset) in 
	lDtObjcts. sTrgtPntr is a set of values for the pointer pntr to 
	extract from lDtObjcts. aprioriMdl is a model selected a priori. 
	aprioriParams are the parameters for the aprioriMdl. If vrbose is 
	True the log message is going to be printed.
	Output: ftdMdlPrmtrs is the frozen model and parameters that 
	maximizes the likelihood of the length data. messg is the log out
	message of the fitting process.
	NOTE: aprioriParams and aprioriMdl shall follow scipy standard 
	names.
	NOTE: mdl is the fitest model out of 81 continous tested by maximum 
	likelihood (by AICc). params are the parameters of the fitest model 
	choose by maximum likelihood. The following information is printed 
	aditionally: lnLklhd is the log likelihood if mdl and params. AICc 
	is the Aikaike information criteria corrected for mdl. D is the D 
	value for a Kolmogorov-Smirnov test for goodness of fit of mdl to 
	the lengths in lDtObjcts, and p is the probability of D. 
	"""
	#----------------------------
	#Assert subset in lDtObjcts
	if sTrgtPntr:
		assert pntr and sTrgtPntr
		try:
			sDiff = sTrgtPntr.difference(set([getattr(dt,pntr) for dt in \
			lDtObjcts]))		
			if sDiff:
				raise exceptions.CelleryWarningObjct \
				('Number of input object not present in lDtObjcts: ', \
				len(sDiff))
		except exceptions.CelleryWarningObjct as err:
				print err
				pass
	#----------------------------
	#Retrieve length for objects with input pointer
	if sTrgtPntr:
		lLenDt = [getattr(dt,'len') for dt in lDtObjcts if \
		getattr(dt,pntr) in sTrgtPntr]
	else:
		lLenDt = [getattr(dt,'len') for dt in lDtObjcts]
	#----------------------------
	#Retrieve model for the input list
	ftdMdlPrmtrs,messg = fitMdltoFlts(lLenDt,aprioriMdl,aprioriParams, \
	outPlt,vrbose)
	return ftdMdlPrmtrs,messg


########################################################
#~ Fit model to a list of float values
########################################################
def fitMdltoFlts(lLenDt,aprioriMdl=False,aprioriParams=False,outPlt=False, \
	vrbose=True):
	"""
	Input: lLenDt is a list with the lengths of the data of interest.
	aprioriMdl is a model selected a priori. aprioriParams are the 
	parameters for the aprioriMdl. outPlt is a file name to plot the 
	distribution of lengths and fitted model. If vrbose is True the log 
	message is going to be printed.
	Output: ftdMdlPrmtrs is the frozen model and parameters that 
	maximizes the likelihood of the length data. messg is the log out
	message of the fitting process.
	NOTE: aprioriParams and aprioriMdl shall follow scipy standard 
	names.
	NOTE: Only values with Kolmogorov-Smirnov p-val>0.05 are considered.
	To avoid overfitting, AICc is the desicion criteria.
	"""
	if aprioriMdl or aprioriParams:
		assert aprioriMdl and aprioriParams
	#----------------------------
	#Define models
	if aprioriMdl:
		aMdlNms = array([aprioriMdl])
	else:
		aMdlNms = array(['alpha','anglit','arcsine','beta','betaprime', \
		'bradford','burr','cauchy','chi','chi2','cosine','dgamma', \
		'dweibull','expon','exponpow','exponweib','f','fatiguelife', \
		'fisk','foldcauchy','foldnorm','frechet_l','frechet_r','gamma', \
		'gausshyper','genexpon','genextreme','gengamma','genhalflogistic', \
		'genlogistic','genpareto','gilbrat','gompertz','gumbel_l', \
		'gumbel_r','halfcauchy','halflogistic','halfnorm','hypsecant', \
		'invgamma','invgauss','invweibull','johnsonsb','johnsonsu', \
		'ksone','kstwobign','laplace','loggamma','logistic','loglaplace', \
		'lognorm','lomax','maxwell','mielke','nakagami','ncf','nct', \
		'ncx2','norm','pareto','pearson3','powerlaw','powerlognorm', \
		'powernorm','rayleigh','rdist','recipinvgauss','reciprocal', \
		'rice','semicircular','t','triang','truncexpon','truncnorm', \
		'tukeylambda','uniform','vonmises','wald','weibull_max', \
		'weibull_min','wrapcauchy'])
	#----------------------------
	#Set required variables
	cntLns = len(lLenDt) # number of input lengths
	nMlds = len(aMdlNms) # number of models
	aMdlAICc = zeros(nMlds,dtype=float32) # holds AICc value
	aMdlD = zeros(nMlds,dtype=float32) # holds D value
	aMdlp = zeros(nMlds,dtype=float32) # holds P-value
	aMdlLnLklhd = zeros(nMlds,dtype=float32) # holds log likelihood
	aMdlfitPrms = zeros(nMlds,dtype=object) # # holds model params.
	aMdlAICc.fill(nan)
	aMdlD.fill(nan)
	aMdlp.fill(nan)
	aMdlLnLklhd.fill(nan)
	aMdlfitPrms.fill(nan)
	#----------------------------
	#Run model fitting for each test and the maximum likelihood model. 
	#Try AICc otherwise try D.
	minAICcVal = inf#set initial minimum abs AICc value
	minAICcPos = nan#set initial minimum abs AICc position
	minDVal = inf#set initial minimum D value
	minDPos = nan#set initial minimum D position
	for mdlPos in xrange(nMlds):
		addVls = True
		aMdlObjct = aMdlNms[mdlPos]
		try:
			mdl = getattr(stats,aMdlObjct)
			AICc,D,p,lnLklhd,params = clcMdl(lLenDt,cntLns,mdl, \
			aprioriParams)
			AICc = absolute(AICc)#absolute AICc value
		except:
			addVls = False
		if addVls:
			aMdlAICc[mdlPos] = AICc
			aMdlD[mdlPos] = D
			aMdlp[mdlPos] = p
			aMdlLnLklhd[mdlPos] = lnLklhd
			aMdlfitPrms[mdlPos] = params
			if p>0.05:
				if AICc<minAICcVal:
					minAICcVal,minAICcPos = AICc,mdlPos
			if D<minDVal:
				minDVal,minDPos = D,mdlPos
	#----------------------------
	#Select the maximum likelihood model. Try AICc otherwise try D.
	messg = []#holds log of the model testing
	if not isnan(minAICcPos):
		AICc = minAICcVal
		D = aMdlD[minAICcPos]
		p = aMdlp[minAICcPos]
		lnLklhd = aMdlLnLklhd[minAICcPos]
		params = aMdlfitPrms[minAICcPos]
		mdl = getattr(stats,aMdlNms[minAICcPos])
	elif not isnan(minDPos):
		AICc = aMdlAICc[minDPos]
		D = minDVal
		p = aMdlp[minDPos]
		lnLklhd = aMdlLnLklhd[minDPos]
		params = aMdlfitPrms[minDPos]
		mdl = getattr(stats,aMdlNms[minDPos])
		messg.extend([
		'\t-----------------------------------------------------------', \
		'\t...Model was choosen but Kolmogorov-Smirnov test %s'% \
		'probability was not significant! The model choose minimizes D', \
		'\t-----------------------------------------------------------'])
	else:
		mdl = getattr(stats,'uniform')
		AICc,D,p,lnLklhd,params = clcMdl(lLenDt,cntLns,mdl, \
		aprioriParams=False)
		try:
			raise exceptions.CelleryWarningObjct \
			('Sample size is too small to estimate parameters:%s,'%
			cntLns,'. A uniform distribution in going to be enforced.')
		except exceptions.CelleryWarningObjct as err:
			print err
			messg.extend([
			'\t-----------------------------------------------------------', \
			'\t%s'%err, \
			'\t-----------------------------------------------------------'])
			pass
	#----------------------------
	#Print message and plot figures
	messg.extend(
	['\t-----------------------------------------------------------',
	'\t   The model choosen to minimize the AICc value was: %s'% mdl.name, \
	'\t   ...with parameters: %s'%','.join([str(p) for p in params]), \
	'\t   ...and AICc value: %s'%AICc, \
	'\t   The ln(likelihood) of the model given the data was: %s'%lnLklhd, \
	'\t   Kolmogorov-Smirnov test D for the data and the %s'% \
	'expected model distribution was: %s'%D, \
	'\t   Kolmogorov-Smirnov test p-value for the data and the %s'% \
	'expected model distribution was: %s'%p, \
	'\t-----------------------------------------------------------'])
	#----------------------------
	#Set output
	ftdMdlPrmtrs = mdl(*params)
	messg = '\n'.join(messg)
	#----------------------------
	#Print message
	if vrbose:#print message
		print messg
	#----------------------------
	#Plot figure
	if outPlt:
		fig = plt.figure()
		nbins = cntLns
		n,bins,patches = plt.hist(lLenDt,histtype='stepfilled', \
		bins=nbins,normed=True,alpha=0.2,label='%s samples'%nbins)
		xLnspc = linspace(min(lLenDt),max(lLenDt),nbins)
		plt.plot(xLnspc,mdl(*params).pdf(xLnspc),label=mdl.name, \
		color='r')
		plt.legend()
		plt.xlabel('Sequence length')
		plt.ylabel('Probability')
		plt.ylim((0,max(n)))
		plt.gcf().subplots_adjust(bottom=0.15)
		plt.savefig(outPlt,bbox_inches='tight',format='eps')
		plt.close()		
	return ftdMdlPrmtrs,messg
	
	
########################################################
#~ Retrieve length models for two objects with not fully empty values in 
# a 2x2 array
########################################################
def rtrnMdlsObjcts2x2NonanVls(aVls,lDtObjctRows,lDtObjctClmns, \
	outPltRows=False,outPltClmns=False,pntr=False,sTrgtPntrRows=False, \
	sTrgtPntrClmns=False):
	"""
	Input: aVls is a 2x2 array with numeric values. lDtObjctRows is a 
	list with gene/lncrna objects with pointers "len" and "pos" for 
	genes/lncrnas of interest with the same positions as rows in aVls. 
	lDtObjctClmns is a list with gene/lncrna objects with pointers "len" 
	and "pos" for genes/lncrnas of interest with the same positions as 
	columns in aVls. outPltRows is a file name to plot the distribution 
	of lengths and fitted model for objects in rows. outPltClmns is a 
	file name to plot the distribution of lengths and fitted model for 
	objects in columns. pntr is the pointer to the attribute of interest 
	(to select a subset) in lDtObjcts. sTrgtPntrRows is a set of values 
	for the pointer pntr to extract from lDtObjctRows and their 
	equivalent positions in aVls. sTrgtPntrClmns is a set of values for 
	the pointer pntr to extract from lDtObjctClmns and their equivalent 
	positions in aVls. 
	Output: ftdMdlPrmtrsRows is the frozen model and parameters that 
	maximizes the likelihood of the length data in rows. messgRows is 
	the log out message of the fitting process for rows. 
	ftdMdlPrmtrsClmns is the frozen model and parameters that maximizes 
	the likelihood of the length data in columns. messgsClmns is the log 
	out message of the fitting process for columns. aNonNullPosRows is 
	an array of positions of non-fully null values in rows. 
	aNonNullPosClmns is an arary of positions of non-fully null values 
	in columns.
	NOTE: if sTrgtPntrRows or sTrgtPntrClmns are true, only non-null 
	values in the subarray made by the two subsets of positions is going
	to be considered.
	"""
	#----------------------------
	#Test for correct size of variables
	lenDtObjctRows,lenDtObjctClmns = aVls.shape
	assert lenDtObjctRows==len(lDtObjctRows)
	assert lenDtObjctClmns==len(lDtObjctClmns)
	if sTrgtPntrRows:
		assert pntr and sTrgtPntrRows
		try:
			sDiff = sTrgtPntrRows.difference(set([getattr(dt,pntr) for \
			dt in lDtObjctRows]))		
			if sDiff:
				raise exceptions.CelleryWarningObjct \
				('Number of input object not present in lDtObjctRows: ', \
				len(sDiff))
		except exceptions.CelleryWarningObjct as err:
				print err
				pass
	if sTrgtPntrClmns:
		assert pntr and sTrgtPntrClmns
		try:
			sDiff = sTrgtPntrClmns.difference(set([getattr(dt,pntr) for \
			dt in lDtObjctClmns]))		
			if sDiff:
				raise exceptions.CelleryWarningObjct \
				('Number of input object not present in lDtObjctClmns: ', \
				len(sDiff))
		except exceptions.CelleryWarningObjct as err:
				print err
				pass
	#----------------------------
	#Retrieve length for objects with input pointer
	aVlsTrgt = copy(aVls)
	if sTrgtPntrRows:
		lPosDtRows = [getattr(dt,'pos') for dt in lDtObjctRows if \
		getattr(dt,pntr) in sTrgtPntrRows]
		#assume sorted(lPosDtRows)==lPosDtRows
		aVlsTrgt = aVlsTrgt[lPosDtRows,:]
	else:
		lPosDtRows = [getattr(dt,'pos') for dt in lDtObjctRows]
	if sTrgtPntr:
		lPosDtClmns = [getattr(dt,'pos') for dt in lDtObjctClmns if \
		getattr(dt,pntr) in sTrgtPntrClmns]
		aVlsTrgt = aVlsTrgt[:,lPosDtClmns]
	else:
		lPosDtClmns = [getattr(dt,'pos') for dt in lDtObjctClmns]
	#----------------------------
	#Retrieve non-null positions in array and models for them
	lPosRowsNonNull,lPosClmnsNonNull = rtrnPos2x2NonanVls(aVlsTrgt)
	ftdMdlPrmtrsRows,messgRows = fitMdl(lDtObjctRows,outPltRows, \
	pntr='pos',sTrgtPntr=set(lPosRowsNonNull))
	ftdMdlPrmtrsClmns,messgsClmns = fitMdl(lDtObjctClmns,outPltClmns, \
	pntr='pos',sTrgtPntr=set(lPosClmnsNonNull))
	aNonNullPosRows = array(lPosRowsNonNull)
	aNonNullPosClmns = array(lPosClmnsNonNull)
	return ftdMdlPrmtrsRows,messgRows,ftdMdlPrmtrsClmns,messgsClmns, \
	aNonNullPosRows,aNonNullPosClmns


########################################################
#~ Retrieve positions for objects with not fully empty values in a 2x2 
# array
########################################################
def rtrnPos2x2NonanVls(aVls):
	"""
	Input: aVls is a 2x2 array with numeric values.
	Output: lPosRowsNonNull is the list of positions in rows of aVls 
	that have at least one non-null value. lPosClmnsNonNull is the list 
	of positions in columns of aVls that have at least one non-null 
	value.
	"""
	#--------------------------------------------------------
	#~ Return non-empty variables and their names
	nonNullMskdAVls = ma.masked_invalid(aVls)
	lenRows,lenClmns = nonNullMskdAVls.shape
	sumRowsNonNullMskdAVls = nonNullMskdAVls.count(1)
	sumClmnsNonNullMskdAVls = nonNullMskdAVls.count(0)
	lPosRowsNonNull = [pos for pos in xrange(lenRows) if \
	sumRowsNonNullMskdAVls[pos]>0]
	lPosClmnsNonNull = [pos for pos in xrange(lenClmns) if \
	sumClmnsNonNullMskdAVls[pos]>0]
	return lPosRowsNonNull,lPosClmnsNonNull



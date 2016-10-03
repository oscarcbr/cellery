#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  randomz.py part of cellery (ceRNAs linking inference)
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
from cellery import exceptions,statistics
from numpy import array,float32,divide,inf,ma,multiply,nan,random, \
vectorize,zeros
from random import seed
from scipy import integrate
from scipy.stats import mstats,binom,norm

import numpy as np


########################################################
#~ Compute randomization on only columns 
########################################################
def cmptClmRndmztn(aBckgrndMetrcMskd,aBckgrndPosRows, \
	aBckgrndClmnPosProb,aTrgtMetrcMskd=None,statstcTrgt=False, \
	stdTrgt=False,lenTrgtRows=False,lenTrgtClms=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosRowsToGns=None,aPosClmnsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None,mirnDtype='cnt', \
	outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. aBckgrndPosRows is the position of 
	the rows in the background. aBckgrndClmnPosProb is the probability 
	of all columns to be selected from the background. Optionally, 
	aTrgtMetrcMskd is a masked array with metric values of interest for 
	the target. statstcTrgt is the statistic value of the target whose 
	probability is going to be calculated from the randomized background 
	using a z-score approach. stdTrgt is the standard deviation of the 
	target data. lenTrgtRows is the number of rows to be sampled. 
	lenTrgtClms is the number of columns to be sampled. outPltFl is a 
	file to plot the randomization and significance of target statistic. 
	numRndmztns is the number of randomizations to run in the background. 
	maxRndmztns is the maximum number of iterations to enforce normality. 
	seedAdd is the seed to run the randomizations. statstc is the 
	statistic to sample from each randomization to build the normal 
	distribution to test the significance of the target statistic. If 
	vrbse all log messages are going to be printed. aPosClmnsToGns is an 
	array of position of gene to which each column in aBckgrndMetrcMskd 
	and aPosProbClmns is mapped. aPosRowsToGns is an array of position 
	of gene to which each row in aBckgrndMetrcMskd and aPosProbRows is 
	mapped. If aPosRowsToGns and aPosClmnsToGns arenot None calculations 
	are going to be run by gene. aBckgrndMirnMetrcMskd is a masked array 
	(3D supermatrix) with an array of metric values for each column and 
	row to be summarized using the value statstc. aBckgrndMirnMetrcMskd 
	is an array with the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names. fnctn is the statistical function to 
	calculate the probability {'sf' is greater or equal than}.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only columns are going to be radomized.
	NOTE: Positions in aBckgrndPosRows are going to be obtained in the 
	same order for every randomization.
	"""
	#----------------------------
	#Test for inputs
	ovrAllogRun = []#holder for log message
	fnctn='sf'#Greater or equal, as normal distributed 1-p is less than.
	assert aTrgtMetrcMskd is not None or (statstcTrgt and stdTrgt and \
	lenTrgtClms)
	assert not lenTrgtRows#no length target for rows
	#mirna inputs
	if aBckgrndMirnMetrcMskd is not None or aTrgtMirnStatstcMskd is not \
		None or outMirnStstFl is not None or aMirnNms is not None:
		assert aBckgrndMirnMetrcMskd is not None and \
		aTrgtMirnStatstcMskd is not None and outMirnStstFl is not None \
		and aMirnNms is not None
		assert mirnDtype in {'cnt','scr'}
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run for','miRNAs')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
	#gene inputs
	if aPosClmnsToGns is not None and aPosRowsToGns is not None:
		mthdBckgrndRndmztn = rnClmRndmztnByGns
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run at the level of','gene')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
		if aTrgtMetrcMskd is not None:
			#calculate stat for the target metrics by gene
			vctrzdRtrnAvrgMskdArray=vectorize(rtrnAvrgMskdArray, \
			excluded={0})
			aTrgtMetrcMskdGns = vctrzdRtrnAvrgMskdArray(aTrgtMetrcMskd, \
			aPosRowsToGns,aPosClmnsToGns)
			aTrgtMetrcMskdGns = ma.masked_invalid(aTrgtMetrcMskdGns)
			statstcTrgt = getattr(ma,statstc)(aTrgtMetrcMskdGns)
			stdTrgt = ma.std(aTrgtMetrcMskdGns)
			lenTrgtClms = aTrgtMetrcMskd.shape[1]
		else:
			assert statstcTrgt and stdTrgt and lenTrgtClms
		if aTrgtMirnStatstcMskd is not None:
			rtrnRndmBckgrndMirnMetrcPrGn = \
			vectorize(rtrnStstMskdMirnArray,excluded={0})
			aTrgtMirnStatstcMskd = \
			rtrnRndmBckgrndMirnMetrcPrGn(aTrgtMirnStatstcMskd, \
			aPosRowsToGns,aPosClmnsToGns,numMrns)
	else:
		mthdBckgrndRndmztn = rnClmRndmztn#only randomize columns
		if statstcTrgt or stdTrgt or lenTrgtClms:
			assert statstcTrgt and stdTrgt and lenTrgtClms
		else:
			statstcTrgt = getattr(ma,statstc)(ma.masked_invalid \
			(aTrgtMetrcMskd))
			stdTrgt = ma.std(ma.masked_invalid(aTrgtMetrcMskd))
			lenTrgtClms = aTrgtMetrcMskd.shape[1]
	#----------------------------
	#Report variables to use
	if seedAdd:
		seedMssg = '\tSeed: %s'%seedAdd
	else:
		seedMssg = ''
	messg= '\n'.join([
	'-----------------------------------------------------------', \
	'A randomization of columns is going to be conducted with parameters:', \
	'\tTarget %s: %s, and standard error: %s'%(statstc,statstcTrgt, \
	stdTrgt), \
	'\t%s columns are going to be sampled in %s randomizations' \
	%(lenTrgtClms,numRndmztns), \
	"\tand significance of target's %s is going to be determined."% \
	statstc, \
	seedMssg, \
	'-----------------------------------------------------------'])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Run randomizations
	p_val,zscore,meanBckgrnd,stdBckgrnd,logRun = rndmzCore \
	(aBckgrndMetrcMskd,mthdBckgrndRndmztn,statstcTrgt,stdTrgt,outPltFl, \
	aBckgrndClmnPosProb,aBckgrndPosRows,lenTrgtClms,lenTrgtRows, \
	seedAdd,numRndmztns,maxRndmztns,statstc,fnctn,vrbse,aPosClmnsToGns, \
	aPosRowsToGns,aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype, \
	outMirnStstFl,aMirnNms)
	#----------------------------
	#Set log out
	if p_val<=0.05:
		signfcn = 'in the top 5% of the values in the randomized'
	elif p_val<=0.01:
		signfcn = 'in the top 1% of the values in the randomized'
	elif p_val>=0.95:
		signfcn = 'in the bottom 5% of the values in the randomized'
		p_val = 1-p_val#p-val in the bottom
	elif p_val>=0.99:
		p_val = 1-p_val#p-val in the bottom
		signfcn = 'in the bottom 1% of the values in the randomized'
	else:
		signfcn = 'to be comparable with the values in the randomized'
	ovrAllogRun.append(logRun)
	messg = '\n'.join([
	'-----------------------------------------------------------', \
	'\t%s randomizations finished!'%numRndmztns, \
	"\tTarget's %s was found %s background"%(statstc,signfcn), \
	"\twith z-score of: %s, and p-value %s."%(zscore,p_val), \
	"\tBackground had normal distribution with mean: %s and std. %s" \
	%(meanBckgrnd,stdBckgrnd)])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Return log message and p-values
	return ovrAllogRun,p_val
	

########################################################
#~ Compute full randomization on an array
########################################################
def cmptFullRndmztn(aBckgrndMetrcMskd,aBckgrndRowsPosProb, \
	aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt,stdTrgt,lenTrgtRows, \
	lenTrgtClms,outPltFl,numRndmztns,maxRndmztns,seedAdd,statstc,vrbse, \
	aPosRowsToGns=None,aPosClmnsToGns=None,aBckgrndMirnMetrcMskd=None, \
	aTrgtMirnStatstcMskd=None,mirnDtype='cnt',outMirnStstFl=None, \
	aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. aBckgrndRowsPosProb is the 
	probability of all rows in the background.  aBckgrndClmnPosProb is 
	the position probability of all columns in the background. 
	Optionally aTrgtMetrcMskd is a masked array with metric values of 
	interest for the target. statstcTrgt is the statistic value of the 
	target whose probability is going to be calculated from the 
	randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled. lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum
	number of iterations to enforce normality. seedAdd is the seed
	to run the randomizations. statstc is the statistic to sample from 
	each randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. If 
	aPosRowsToGns and aPosClmnsToGns are not None calculations are going 
	to be run by gene. aBckgrndMirnMetrcMskd is a masked array (3D 
	supermatrix) with an array of metric values for each column and row 
	to be summarized using the value statstc. aBckgrndMirnMetrcMskd is 
	an array with the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Both, rows and columns are going to be radomized.
	"""
	#----------------------------
	#Test for inputs
	ovrAllogRun = []#holder for log message
	fnctn='sf'#Greater or equal, as normal distributed 1-p is less than.
	assert aTrgtMetrcMskd is not None or (statstcTrgt and stdTrgt and \
	lenTrgtRows and lenTrgtClms)
	#mirna inputs
	if aBckgrndMirnMetrcMskd is not None or aTrgtMirnStatstcMskd is not \
		None or outMirnStstFl is not None or aMirnNms is not None:
		assert aBckgrndMirnMetrcMskd is not None and \
		aTrgtMirnStatstcMskd is not None and outMirnStstFl is not None \
		and aMirnNms is not None
		assert mirnDtype in {'cnt','scr'}
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run for','miRNAs')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
	#gene inputs
	if aPosClmnsToGns is not None and aPosRowsToGns is not None:
		mthdBckgrndRndmztn = rnClmNRowRndmztnByGns
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run at the level of','gene')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
		if aTrgtMetrcMskd is not None:
			#calculate stat for the target metrics by gene
			vctrzdRtrnAvrgMskdArray = vectorize(rtrnAvrgMskdArray, \
			excluded={0})
			aTrgtMetrcMskdGns = vctrzdRtrnAvrgMskdArray(aTrgtMetrcMskd, \
			aPosRowsToGns,aPosClmnsToGns)
			aTrgtMetrcMskdGns = ma.masked_invalid(aTrgtMetrcMskdGns)
			statstcTrgt = getattr(ma,statstc)(aTrgtMetrcMskdGns)
			stdTrgt = ma.std(aTrgtMetrcMskdGns)
			lenTrgtRows,lenTrgtClms = aTrgtMetrcMskd.shape
		else:
			assert statstcTrgt and stdTrgt and lenTrgtRows and \
			lenTrgtClms
		if aTrgtMirnStatstcMskd is not None:
			rtrnRndmBckgrndMirnMetrcPrGn = \
			vectorize(rtrnStstMskdMirnArray,excluded={0})
			aTrgtMirnStatstcMskd = \
			rtrnRndmBckgrndMirnMetrcPrGn(aTrgtMirnStatstcMskd, \
			aPosRowsToGns,aPosClmnsToGns,numMrns)
	else:
		mthdBckgrndRndmztn = rnClmNRowRndmztn#full randomz. method.
		if statstcTrgt or stdTrgt or lenTrgtRows or lenTrgtClms:
			assert statstcTrgt and stdTrgt and lenTrgtRows and \
			lenTrgtClms
		else:
			statstcTrgt = getattr(ma,statstc)(ma.masked_invalid \
			(aTrgtMetrcMskd))
			stdTrgt = ma.std(ma.masked_invalid(aTrgtMetrcMskd))
			lenTrgtRows,lenTrgtClms = aTrgtMetrcMskd.shape
	#----------------------------
	#Report variables to use
	if seedAdd:
		seedMssg = '\tSeed: %s'%seedAdd
	else:
		seedMssg = ''
	messg= '\n'.join([
	'-----------------------------------------------------------', \
	'A full randomization is going to be conducted with parameters:', \
	'\tTarget %s: %s, and standard error: %s'%(statstc,statstcTrgt, \
	stdTrgt), \
	'\t%s columns and %s rows are going to be sampled in %s randomizations' \
	%(lenTrgtClms,lenTrgtRows,numRndmztns), \
	"\tand significance of target's %s is going to be determined."% \
	statstc, \
	seedMssg, \
	'-----------------------------------------------------------'])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Run randomizations
	p_val,zscore,meanBckgrnd,stdBckgrnd,logRun = rndmzCore \
	(aBckgrndMetrcMskd,mthdBckgrndRndmztn,statstcTrgt,stdTrgt,outPltFl, \
	aBckgrndClmnPosProb,aBckgrndRowsPosProb,lenTrgtClms,lenTrgtRows, \
	seedAdd,numRndmztns,maxRndmztns,statstc,fnctn,vrbse,aPosClmnsToGns, \
	aPosRowsToGns,aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype, \
	outMirnStstFl,aMirnNms)
	#----------------------------
	#Set log out
	if p_val<=0.05:
		signfcn = 'in the top 5% of the values in the randomized'
	elif p_val<=0.01:
		signfcn = 'in the top 1% of the values in the randomized'
	elif p_val>=0.95:
		signfcn = 'in the bottom 5% of the values in the randomized'
		p_val = 1-p_val#p-val in the bottom
	elif p_val>=0.99:
		p_val = 1-p_val#p-val in the bottom
		signfcn = 'in the bottom 1% of the values in the randomized'
	else:
		signfcn = 'to be comparable with the values in the randomized'
	ovrAllogRun.append(logRun)
	messg = '\n'.join([
	'-----------------------------------------------------------', \
	'\t%s randomizations finished!'%numRndmztns, \
	"\tTarget's %s was found %s background"%(statstc,signfcn), \
	"\twith z-score of: %s, and p-value %s."%(zscore,p_val), \
	"\tBackground had normal distribution with mean: %s and std. %s" \
	%(meanBckgrnd,stdBckgrnd)])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Return log message and p-values
	return ovrAllogRun,p_val


########################################################
#~ Compute randomization on only rows 
########################################################
def cmptRowRndmztn(aBckgrndMetrcMskd,aBckgrndRowsPosProb, \
	aBckgrndPosClmns,aTrgtMetrcMskd=None,statstcTrgt=False,stdTrgt=False, \
	lenTrgtRows=False,lenTrgtClms=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosRowsToGns=None,aPosClmnsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None,mirnDtype='cnt', \
	outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. aBckgrndRowsPosProb is the 
	probability of all rows in the background. aBckgrndPosClmns is the 
	position of the columns to be selected from the background. 
	Optionally, aTrgtMetrcMskd is a masked array with metric values of 
	interest for the target. statstcTrgt is the statistic value of the 
	target whose probability is going to be calculated from the 
	randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled. lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum
	number of iterations to enforce normality. seedAdd is the seed
	to run the randomizations. statstc is the statistic to sample from 
	each randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. If 
	aPosRowsToGns and aPosClmnsToGns are not None calculations are going 
	to be run by gene. aBckgrndMirnMetrcMskd is a masked array (3D 
	supermatrix) with an array of metric values for each column and row 
	to be summarized using the value statstc. aBckgrndMirnMetrcMskd is 
	an array with the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only rows are going to be radomized.
	NOTE: Positions in aBckgrndPosClmns are going to be obtained in the 
	same order for every randomization.
	"""
	#----------------------------
	#Test for inputs
	ovrAllogRun = []#holder for log message
	fnctn='sf'#Greater or equal, as normal distributed 1-p is less than.
	assert aTrgtMetrcMskd is not None or (statstcTrgt and stdTrgt and \
	lenTrgtRows)
	assert not lenTrgtClms#no length target for rows
	#mirna inputs
	if aBckgrndMirnMetrcMskd is not None or aTrgtMirnStatstcMskd is not \
		None or outMirnStstFl is not None or aMirnNms is not None:
		assert aBckgrndMirnMetrcMskd is not None and \
		aTrgtMirnStatstcMskd is not None and outMirnStstFl is not None \
		and aMirnNms is not None
		assert mirnDtype in {'cnt','scr'}
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run for','miRNAs')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
	#gene inputs
	if aPosClmnsToGns is not None and aPosRowsToGns is not None:
		mthdBckgrndRndmztn = rnRowRndmztnByGns
		try:
			raise exceptions.CelleryWarningObjct \
			('Randomizations are going to be run at the level of','gene')
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			ovrAllogRun.append(mssge)
			pass
		if aTrgtMetrcMskd is not None:
			#calculate stat for the target metrics by gene
			vctrzdRtrnAvrgMskdArray=vectorize(rtrnAvrgMskdArray, \
			excluded={0})
			aTrgtMetrcMskdGns = vctrzdRtrnAvrgMskdArray(aTrgtMetrcMskd, \
			aPosRowsToGns,aPosClmnsToGns)
			aTrgtMetrcMskdGns = ma.masked_invalid(aTrgtMetrcMskdGns)
			statstcTrgt = getattr(ma,statstc)(aTrgtMetrcMskdGns)
			stdTrgt = ma.std(aTrgtMetrcMskdGns)
			lenTrgtRows = aTrgtMetrcMskd.shape[0]
		else:
			assert statstcTrgt and stdTrgt and lenTrgtRows
		if aTrgtMirnStatstcMskd is not None:
			rtrnRndmBckgrndMirnMetrcPrGn = \
			vectorize(rtrnStstMskdMirnArray,excluded={0})
			aTrgtMirnStatstcMskd = \
			rtrnRndmBckgrndMirnMetrcPrGn(aTrgtMirnStatstcMskd, \
			aPosRowsToGns,aPosClmnsToGns,numMrns)
	else:
		mthdBckgrndRndmztn = rnRowRndmztn#only randomize rows
		if statstcTrgt or stdTrgt or lenTrgtRows:
			assert statstcTrgt and stdTrgt and lenTrgtRows
		else:
			statstcTrgt = getattr(ma,statstc)(ma.masked_invalid \
			(aTrgtMetrcMskd))
			stdTrgt = ma.std(ma.masked_invalid(aTrgtMetrcMskd))
			lenTrgtRows = aTrgtMetrcMskd.shape[0]
	#----------------------------
	#Report variables to use
	if seedAdd:
		seedMssg = '\tSeed: %s'%seedAdd
	else:
		seedMssg = ''
	messg= '\n'.join([
	'-----------------------------------------------------------', \
	'A randomization of rows is going to be conducted with parameters:', \
	'\tTarget %s: %s, and standard error: %s'%(statstc,statstcTrgt, \
	stdTrgt), \
	'\t%s rows are going to be sampled in %s randomizations' \
	%(lenTrgtRows,numRndmztns), \
	"\tand significance of target's %s is going to be determined."% \
	statstc, \
	seedMssg, \
	'-----------------------------------------------------------'])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Run randomizations
	p_val,zscore,meanBckgrnd,stdBckgrnd,logRun = rndmzCore \
	(aBckgrndMetrcMskd,mthdBckgrndRndmztn,statstcTrgt,stdTrgt,outPltFl, \
	aBckgrndPosClmns,aBckgrndRowsPosProb,lenTrgtClms,lenTrgtRows, \
	seedAdd,numRndmztns,maxRndmztns,statstc,fnctn,vrbse,aPosClmnsToGns, \
	aPosRowsToGns,aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype, \
	outMirnStstFl,aMirnNms)
	#----------------------------
	#Set log out
	if p_val<=0.05:
		signfcn = 'in the top 5% of the values in the randomized'
	elif p_val<=0.01:
		signfcn = 'in the top 1% of the values in the randomized'
	elif p_val>=0.95:
		signfcn = 'in the bottom 5% of the values in the randomized'
		p_val = 1-p_val#p-val in the bottom
	elif p_val>=0.99:
		p_val = 1-p_val#p-val in the bottom
		signfcn = 'in the bottom 1% of the values in the randomized'
	else:
		signfcn = 'to be comparable with the values in the randomized'
	ovrAllogRun.append(logRun)
	messg = '\n'.join([
	'-----------------------------------------------------------', \
	'\t%s randomizations finished!'%numRndmztns, \
	"\tTarget's %s was found %s background"%(statstc,signfcn), \
	"\twith z-score of: %s, and p-value %s."%(zscore,p_val), \
	"\tBackground had normal distribution with mean: %s and std. %s" \
	%(meanBckgrnd,stdBckgrnd)])
	ovrAllogRun.append(messg)
	if vrbse:
		print messg
	#----------------------------
	#Return log message and p-values
	return ovrAllogRun,p_val
	

########################################################
#~ Calculate a statistic on miRNA metrics for a pair of genes/lncrnas.
########################################################
def cmptMirnaStatRowClmnPrs(aRowsMirnMetrc,aClmnsMirnMetrc, \
	aRowsClmnsMsk=None,mirnPrsStstc='add',dtype='cnt'):
	"""
	Input: aRowsMirnMetrc is an array of miRNA metric values arrays for 
	each row. aClmnsMirnMetrc is an array of miRNA metric values arrays 
	for each column. Optionally, aRowsClmnsMsk is an array of size 
	len(aRowsMirnMetrc) x len(aClmnsMirnMetrc)to mask the output matrix.
	mirnPrsStstc is the statistic to calculate values for miRNAs. dtype
	is the datatype of the miRNA metric: {'cnt' for counts and 'scr' 
	for scores}.
	Output: aMirnMetrcMskd is an output array with stats for the metrics 
	of shared miRNAs for each masked row-column pair.
	NOTE: row-column pairs masked in aRowsClmnsMsk are retrieved but 
	they are masked.
	NOTE: non-shared miRNAs have values of nan.
	"""
	#----------------------------
	#Test input variables
	lenRows,lenMirnsRows = aRowsMirnMetrc.shape
	lenClmns,lenMirnsClmns = aClmnsMirnMetrc.shape
	assert lenMirnsRows==lenMirnsClmns
	assert dtype in {'cnt','scr'}
	if aRowsClmnsMsk is not None:
		lenRowsMsk,lenClmnsMsk,lenMirnsMsk = aRowsClmnsMsk.shape
		assert lenRows==lenRowsMsk and lenClmnsMsk==lenClmns and \
		lenMirnsMsk==lenMirnsRows
	else:#make an empty mask
		aRowsClmnsMsk = zeros((lenRows,lenClmns,lenMirnsRows))
	if dtype=='scr':
		aRowsMirnMetrcMskd = ma.masked_invalid(aRowsMirnMetrc)
		aClmnsMirnMetrcMskd = ma.masked_invalid(aClmnsMirnMetrc)
	else:
		aRowsMirnMetrcMskd = ma.masked_less(aRowsMirnMetrc,1)
		aClmnsMirnMetrcMskd = ma.masked_less(aClmnsMirnMetrc,1)
	#----------------------------
	#Calculate metric for each miRNA in each row-column pair
	aMirnMetrc = ma.zeros((lenRows,lenClmns,lenMirnsRows),dtype=float32)
	for rowPos in xrange(lenRows):
		for clmnPos in xrange(lenClmns):
			aMirnMetrc[rowPos][clmnPos] = getattr(ma,mirnPrsStstc) \
			(aRowsMirnMetrcMskd[rowPos],aClmnsMirnMetrcMskd[clmnPos])
	aMirnMetrc.fill_value = nan
	#----------------------------
	#Mask the output stats array
	aMirnMetrcMskd = ma.array(aMirnMetrc.filled(),mask=aRowsClmnsMsk)
	return aMirnMetrcMskd


########################################################
#~ Plot data distribution given a background
########################################################
def mkPlt(statstcTrgt,aBcrkngdStatstc,outPltFl):
	"""
	Input: statstcTrgt is the target statistic value. aBcrkngdStatstc is
	an array of statistic values for the background sampling. outPltFl 
	is the output plot file.
	Output: outPltFl is the output plot file with the distribution of 
	the background indicating the position of the target.
	"""
	fig = plt.figure()
	xVls,yVls,info = plt.hist(aBcrkngdStatstc,100,normed=True)
	plt.xlabel('Mean ECP value')
	plt.ylabel('Number of replicates')
	plt.annotate('Target (%s)'%statstcTrgt, xy=(statstcTrgt, max(xVls)), \
	xytext = (statstcTrgt, max(xVls)+10),arrowprops=dict(facecolor='red', \
	shrink = 0.05, width=2))
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.savefig(outPltFl,bbox_inches='tight',format='eps')
	plt.close()	
	return 0	


########################################################
#~ Process miRNA randomization results
########################################################
def procMirnRndmztn(aBcrkngdMirnStatstc,aTrgtMirnStatstcMskd,mirnDtype, \
	outMirnStstFl,aMirnNms,fnctn='sf'):
	"""
	Input: aBcrkngdMirnStatstc is an array with the miRNA metrics 
	summarized for the background. aTrgtMirnStatstcMskd is an array with 
	the miRNA metrics for the target dataset. mirnDtype is the datatype 
	of the miRNA metric: {'cnt' for counts and 'scr' for scores}. 
	outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names. fnctn
	is the statistical function to calculate the probability {'sf' is 
	greater or equal than}.
	Output: mssge is the log text of the probability calculation 
	results.
	NOTE: Results are going to be written in the outMirnStstFl file.
	NOTE: aBcrkngdMirnStatstc is an unmasked array and 
	aTrgtMirnStatstcMskd is masked.
	"""
	#----------------------------
	#Calculate statistic for target dataset and set initial parameters
	aTrgtMirnStatstcMskdFlld = rtrnMirnStat(aTrgtMirnStatstcMskd)	
	aTrgtMirnMskdAvrg = rtrnMirnStat(aTrgtMirnStatstcMskd,'mean')
	aTrgtMirnMskdStd = rtrnMirnStat(aTrgtMirnStatstcMskd,'std')
	aBcrkngdMirnStatstcMskd = ma.masked_invalid(aBcrkngdMirnStatstc)
	aBcrkngdMirnStatstcMskdAvrg = rtrnMirnStat(aBcrkngdMirnStatstcMskd, \
	'mean')
	aBcrkngdMirnStatstcMskdStd = rtrnMirnStat(aBcrkngdMirnStatstcMskd, \
	'std')
	nReps,nMirns = aBcrkngdMirnStatstcMskd.shape
	assert nMirns==len(aTrgtMirnStatstcMskdFlld)==len(aMirnNms)
	#----------------------------
	#Calculate significance of miRNA scores/counts using z-score
	aMirPval = zeros(nMirns,dtype=float32)
	aMirZscore = zeros(nMirns,dtype=float32)
	aMirPval.fill(nan)
	aMirZscore.fill(nan)
	for mirnPos in xrange(nMirns):
		mirTrgtStat = aTrgtMirnStatstcMskdFlld[mirnPos]
		mirBckgrndStat = aBcrkngdMirnStatstcMskd[:,mirnPos]
		mirZscore = mstats.zmap(mirTrgtStat,mirBckgrndStat)
		aMirZscore[mirnPos] = mirZscore
		#equal or greater than if 'sf'
		aMirPval[mirnPos] = getattr(norm,fnctn)(mirZscore)
	aMirQvalGrtrEql = statistics.bh_qvalues(aMirPval)#FDR values gtrEql 
	aMirQvalLssEql = statistics.bh_qvalues(1-aMirPval)#FDR values lessEql 
	#----------------------------
	#Calculate significance of miRNA counts using exact binomial
	if mirnDtype=='cnt':
		#number of trials
		N = ma.sum(aTrgtMirnStatstcMskd).__float__()
		N = float32(N)
		#probability of sucess single trial
		aMirnP_sngl = ma.sum(aBcrkngdMirnStatstcMskd,axis=0)/ma.sum \
		(aBcrkngdMirnStatstcMskd)
		aMirnP_sngl = aMirnP_sngl.filled()
		#equal or greater than
		aMirnCntPvalGrtrEql = binom.sf(aTrgtMirnStatstcMskdFlld-1,N, \
		aMirnP_sngl)
		#equal or less than
		aMirnCntPvalLssEql = 1-binom.sf(aTrgtMirnStatstcMskdFlld,N, \
		aMirnP_sngl)
		aMirnCntQvalGrtrEql = statistics.bh_qvalues(aMirnCntPvalGrtrEql)
		aMirnCntQvalLssEql = statistics.bh_qvalues(aMirnCntPvalLssEql)	
	#----------------------------
	#Format output
	cntSngfMoreZscore = 0
	cntSngfLssZscore = 0
	if mirnDtype=='cnt':
		cntSngfLssBnml = 0
		cntSngfMoreBnml = 0
		outRslts = ['\t'.join(['name','target_average','target_std', \
		'target_statistic','background_statistic_average', \
		'background_statistic_std','target_statistic_k2', \
		'target_statistic_k-wPval','z-score', \
		'p-val (target_statistic >= background_statistic)', \
		'q-val (target_statistic >= background_statistic)', \
		'p-val (target_statistic =< background_statistic)', \
		'q-val (target_statistic =< background_statistic)', \
		'significance (target_statistic >= background_statistic)', \
		'significance (target_statistic =< background_statistic)', \
		'N (binomial test parameter)','k (observations)','p sucess', \
		'p-val (target_counts >= background_counts)', \
		'q-val (target_counts >= background_counts)', \
		'p-val (target_counts =< background_counts)', \
		'q-val (target_counts =< background_counts)', \
		'significance (target_counts >= background_counts)', \
		'significance (target_counts =< background_counts)'])]
	else:
		outRslts = ['\t'.join(['name','target_average','target_std', \
		'target_statistic','background_statistic_average', \
		'background_statistic_std','target_statistic_k2', \
		'target_statistic_k-wPval','z-score', \
		'p-val (target_statistic >= background_statistic)', \
		'q-val (target_statistic >= background_statistic)', \
		'p-val (target_statistic =< background_statistic)', \
		'q-val (target_statistic =< background_statistic)', \
		'significance (target_statistic >= background_statistic)', \
		'significance (target_statistic =< background_statistic)'])]
	for mirnPos in xrange(nMirns):
		mirnNm = aMirnNms[mirnPos]
		mirTrgtArvg = aTrgtMirnMskdAvrg[mirnPos]
		mirTrgtStd = aTrgtMirnMskdStd[mirnPos]
		mirTrgtStat = aTrgtMirnStatstcMskdFlld[mirnPos]
		mirBckgrndStat = aBcrkngdMirnStatstcMskd[:,mirnPos]
		mirBckgrndStatArvg = aBcrkngdMirnStatstcMskdAvrg[mirnPos]
		mirBckgrndStatStd = aBcrkngdMirnStatstcMskdStd[mirnPos]
		mirK2,aMirKWPval = mstats.normaltest(mirBckgrndStat)
		mirZscore = aMirZscore[mirnPos] 
		mirPvalGrtrEql = aMirPval[mirnPos] 
		mirQvalGrtrEql = aMirQvalGrtrEql[mirnPos]
		signfcGrtrEql = ''
		if mirQvalGrtrEql<=0.05:
			signfcGrtrEql+='*'
			cntSngfMoreZscore+=1
		if mirQvalGrtrEql<=0.01:
			signfcGrtrEql+='*'
		mirPvalLssEql = 1-mirPvalGrtrEql
		mirQvalLssEql = aMirQvalLssEql[mirnPos]
		signfcLssEql = ''
		if mirQvalLssEql<=0.05:
			signfcLssEql+='*'
			cntSngfLssZscore+=1
		if mirQvalLssEql<=0.01:
			signfcLssEql+='*'
		if mirnDtype=='cnt':
			signfcCntGrtrEql = ''
			signfcCntLssEql = ''
			mirN = N
			mirK  = aTrgtMirnStatstcMskdFlld[mirnPos]
			mirPsuccs = aMirnP_sngl[mirnPos]
			mirPCntGrtrEql = aMirnCntPvalGrtrEql[mirnPos]
			mirQCntGrtrEql = aMirnCntQvalGrtrEql[mirnPos]
			if mirQCntGrtrEql<=0.05:
				signfcCntGrtrEql+='*'
				cntSngfMoreBnml+=1
			if mirQCntGrtrEql<=0.01:
				signfcCntGrtrEql+='*'
			mirPCntLssEql = aMirnCntPvalLssEql[mirnPos]
			mirQCntLssEql = aMirnCntQvalLssEql[mirnPos]
			if mirQCntLssEql<=0.05:
				signfcCntLssEql+='*'
				cntSngfLssBnml+=1
			if mirQCntLssEql<=0.01:
				signfcCntLssEql+='*'
			outRslts.append('\t'.join([str(v) for v in mirnNm, \
			mirTrgtArvg,mirTrgtStd,mirTrgtStat,mirBckgrndStatArvg, \
			mirBckgrndStatStd,mirK2,aMirKWPval,mirZscore,mirPvalGrtrEql, \
			mirQvalGrtrEql,	mirPvalLssEql,mirQvalLssEql,signfcGrtrEql, \
			signfcLssEql,mirN,mirK,mirPsuccs,mirPCntGrtrEql, \
			mirQCntGrtrEql,mirPCntLssEql,mirQCntLssEql,signfcCntGrtrEql, \
			signfcCntLssEql]))
		else:
			outRslts.append('\t'.join([str(v) for v in mirnNm, \
			mirTrgtArvg,mirTrgtStd,mirTrgtStat,mirBckgrndStatArvg, \
			mirBckgrndStatStd,mirK2,aMirKWPval,mirZscore,mirPvalGrtrEql, \
			mirQvalGrtrEql,	mirPvalLssEql,mirQvalLssEql,signfcGrtrEql, \
			signfcLssEql]))
	#----------------------------
	#Write output file
	outFl = open(outMirnStstFl,'w')
	outFl.write('\n'.join(outRslts))
	outFl.close()
	#----------------------------
	#Make output message
	if mirnDtype == 'cnt':
		mthds = 'z-score and binomial probability'
		sgnfcntRprt = \
		'%s miRNAs for z-score and %s for binomial probability'% \
		(cntSngfMoreZscore,cntSngfMoreBnml)
	else:
		mthds = 'z-score'
		sgnfcntRprt = '%s miRNAs for z-score'%cntSngfMoreZscore
	mssge = '\n'.join([
	'\t-----------------------------------------------------------',
	'\t   The probabilities for independent miRNAs were calculated',
	'\t   using: %s.'%mthds,
	'\t   Results were written in the file: %s'%outMirnStstFl,
	'\t   %s were found to be significant (q<0.05)'%sgnfcntRprt,
	'\t-----------------------------------------------------------'])
	return mssge


########################################################
#~ Run randomization on (both) rows and columns
########################################################
def rnClmNRowRndmztn(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosProbClmns,aPosProbRows,lenTrgtClms=False,lenTrgtRows=False, \
	aPosClmnsToGns=None,aPosRowsToGns=None,aBckgrndMirnMetrcMskd=None, \
	vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed. aPosProbClmns is an array with the 
	probabilities of sampling the positions in columns. aPosProbRows is 
	an array with the probabilities of sampling the positions in rows. 
	lenTrgtClms is the size from the column samples to take (normarlly 
	the same as the target).lenTrgtRows is the size of the samples to 
	take from the rows (normarlly the same as the target). 
	aPosClmnsToGns and aPosRowsToGns are always None. Optionally, 
	aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) with an 
	array of metric values for each column and row to be summarized 
	using the value statstc. 
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of columns and rows (of size lenTrgtClms and 
	lenTrgtRows respectively) in numRndmztns. aBcrkngdMirnStatstc is an 
	with the miRNA metrics summarized, is None if aBckgrndMirnMetrcMskd 
	isNone.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbClmns = len(aPosProbClmns)
	lenAPosProbRows = len(aPosProbRows)
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nRows==lenAPosProbRows and nClmns==lenAPosProbClmns
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		aPosClmns = xrange(lenAPosProbClmns)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosClmsRndm = random.choice(aPosClmns,lenTrgtClms, \
		p=aPosProbClmns)
		#
		aPosRows = xrange(lenAPosProbRows)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosRowsRndm = random.choice(aPosRows,lenTrgtRows, \
		p=aPosProbRows)
		aRndmBckgrndMetrc = aBckgrndMetrcMskd[:,aPosClmsRndm]
		aRndmBckgrndMetrc = aRndmBckgrndMetrc[aPosRowsRndm,:]
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = aBckgrndMirnMetrcMskd[:,aPosClmsRndm]
			aRndmBckgrndMirnMtrc = aRndmBckgrndMirnMtrc[aPosRowsRndm,:]
			aRndmBckgrndMirnMtrc.fill_value = nan
			aRndmBckgrndMirnMtrc = ma.masked_invalid \
			(aRndmBckgrndMirnMtrc.filled())
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Run randomization on (both) rows and columns and average by gene
########################################################
def rnClmNRowRndmztnByGns(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosProbClmns,aPosProbRows,lenTrgtClms,lenTrgtRows, \
	aPosClmnsToGns,aPosRowsToGns,aBckgrndMirnMetrcMskd=None,vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed. aPosProbClmns is an array with the 
	probabilities of sampling the positions in columns. aPosProbRows is 
	an array with the probabilities of sampling the positions in rows. 
	lenTrgtClms is the size from the column samples to take (normarlly 
	the same as the target).lenTrgtRows is the size of the samples to 
	take from the rows (normarlly the same as the target).aPosClmnsToGns 
	is an array of position of gene to which each column in 
	aBckgrndMetrcMskd and aPosProbClmns is mapped. aPosRowsToGns is an 
	array of position of gene to which each row in aBckgrndMetrcMskd and 
	aPosProbRows is mapped. Optionally, aBckgrndMirnMetrcMskd is a 
	masked array (3D supermatrix) with an array of metric values for 
	each column and row to be summarized using the value statstc.
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of columns and rows (of size lenTrgtClms and 
	lenTrgtRows respectively) in numRndmztns. aBcrkngdMirnStatstc is an 
	with the miRNA metrics summarized, is None if aBckgrndMirnMetrcMskd 
	is None.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbClmns = len(aPosProbClmns)
	lenAPosProbRows = len(aPosProbRows)
	#vectorize the function
	rtrnRndmBckgrndMetrcPrGn =vectorize(rtrnAvrgMskdArray,excluded={0})
	lenClmnsToGns = max(aPosClmnsToGns)
	lenRowsToGns = max(aPosRowsToGns)
	#to randomize miRNA metrics
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nRows==lenAPosProbRows and nClmns==lenAPosProbClmns
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
		rtrnRndmBckgrndMirnMetrcPrGn = \
		vectorize(rtrnStstMskdMirnArray,excluded={0})
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		aPosClmns = xrange(lenAPosProbClmns)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosClmsRndmGns = zeros(lenClmnsToGns,dytpe=float32)
		aPosClmsRndmGns.fill(nan)
		aPosRowsRndmGns = zeros(lenRowsToGns,dytpe=float32)
		aPosRowsRndmGns.fill(nan)
		#~ 
		aPosClmsRndm = random.choice(aPosClmns,lenTrgtClms, \
		p=aPosProbClmns)
		for rltvPos,absPos in enumerate(aPosClmsRndm):
			aPosClmsRndmGns[aPosClmnsToGns[absPos]].append(rltvPos)
		#
		aPosRows = xrange(lenAPosProbRows)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosRowsRndm = random.choice(aPosRows,lenTrgtRows, \
		p=aPosProbRows)
		for rltvPos,absPos in enumerate(aPosRowsRndm):
			aPosRowsRndmGns[aPosRowsToGns[absPos]].append(rltvPos)
		#
		aRndmBckgrndMetrc = rtrnRndmBckgrndMetrcPrGn(aBckgrndMetrcMskd, \
		aPosRowsRndmGns,aPosClmsRndmGns)
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = rtrnRndmBckgrndMirnMetrcPrGn \
			(aBckgrndMirnMetrcMskd,aPosRowsRndmGns,aPosClmsRndmGns, \
			numMrns)
			aRndmBckgrndMirnMtrc=ma.masked_invalid(aRndmBckgrndMirnMtrc)
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Run randomization on columns
########################################################
def rnClmRndmztn(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosProbClmns=None,aPosRows=None,lenTrgtClms=False, \
	lenTrgtRows=False,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed. aPosRows is an array of positions in rows to be 
	sampled. aPosProbClmns is an array with the probabilities of 
	sampling the positions in columns. lenTrgtClms is the size of the 
	samples to take from the columns (normarlly the same as the target). 
	aPosClmnsToGns and aPosRowsToGns are always None. Optionally, 
	aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) with an 
	array of metric values for each column and row to be summarized 
	using the value statstc. 
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of columns (of size lenTrgtClms) and values in 
	aPosRows for numRndmztns randomizations. aBcrkngdMirnStatstc is an 
	with the miRNA metrics summarized, is None if aBckgrndMirnMetrcMskd 
	is None.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbClmns = len(aPosProbClmns)
	if aPosRow is None:
		aPosRows = xrange(aBckgrndMetrcMskd.shape[0])
		try:
			raise exceptions.CelleryWarningObjct \
			('Columns are going to be randomized but no positions were',
			'provided for rows, all of them are being included.')
		except exceptions.CelleryWarningObjct as err:
			print err
			pass
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nClmns==lenAPosProbClmns
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		#
		aPosClmns = xrange(lenAPosProbClmns)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosClmsRndm = random.choice(aPosClmns,lenTrgtClms, \
		p=aPosProbClmns)
		aRndmBckgrndMetrc = aBckgrndMetrcMskd[:,aPosClmsRndm]
		aRndmBckgrndMetrc = aRndmBckgrndMetrc[aPosRows,:]
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = aBckgrndMirnMetrcMskd[:,aPosClmsRndm]
			aRndmBckgrndMirnMtrc = aRndmBckgrndMirnMtrc[aPosRows,:]
			aRndmBckgrndMirnMtrc.fill_value = nan
			aRndmBckgrndMirnMtrc = ma.masked_invalid \
			(aRndmBckgrndMirnMtrc.filled())
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Run randomization on only columns and average by gene
########################################################
def rnClmRndmztnByGns(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosProbClmns=None,aPosRows=None,lenTrgtClms=None, \
	lenTrgtRows=None,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed. aPosProbClmns is an array with the 
	probabilities of sampling the positions in columns. aPosClmns is an 
	array ofpositions in rows to be sampled. lenTrgtClms is the size 
	from the column samples to take (normarlly the same as the target). 
	lenTrgtRows is the size of the samples to take from the rows 
	(normally the same as the target).aPosClmnsToGns is an array of 
	position of gene to which each column in aBckgrndMetrcMskd and 
	aPosProbClmns is mapped. aPosRowsToGns is an array of position of 
	gene to which each row in aBckgrndMetrcMskd and aPosRows is mapped.
	Optionally, aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) 
	with an array of metric values for each column and row to be 
	summarized using the value statstc.
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of columns and rows (of size lenTrgtClms and 
	lenTrgtRows respectively) in numRndmztns. aBcrkngdMirnStatstc is an 
	with the miRNA metrics summarized, is None if aBckgrndMirnMetrcMskd 
	is None.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbClmns = len(aPosProbClmns)
	lenClmnsToGns = max(aPosClmnsToGns)
	lenRowsToGns = max(aPosRowsToGns)
	#vectorize the function
	rtrnRndmBckgrndMetrcPrGn =vectorize(rtrnAvrgMskdArray,excluded={0})
	if aPosRow is None:
		aPosRows = xrange(aBckgrndMetrcMskd.shape[0])
		try:
			raise exceptions.CelleryWarningObjct \
			('Columns are going to be randomized but no positions were',
			'provided for rows, all of them are being included.')
		except exceptions.CelleryWarningObjct as err:
			print err
			pass
	#to randomize miRNA metrics
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nClmns==lenAPosProbClmns
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
		rtrnRndmBckgrndMirnMetrcPrGn = \
		vectorize(rtrnStstMskdMirnArray,excluded={0})
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		#
		aPosClmsRndmGns = zeros(lenClmnsToGns,dytpe=float32)
		aPosClmsRndmGns.fill(nan)
		aPosRowsRndmGns = zeros(lenRowsToGns,dytpe=float32)
		aPosRowsRndmGns.fill(nan)
		#~ 
		aPosClmns = xrange(lenAPosProbClmns)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosClmsRndm = random.choice(aPosClmns,lenTrgtClms, \
		p=aPosProbClmns)
		for rltvPos,absPos in enumerate(aPosClmsRndm):
			aPosClmsRndmGns[aPosClmnsToGns[absPos]].append(rltvPos)
		for rltvPos,absPos in enumerate(aPosRows):
			aPosRowsRndmGns[aPosRowsToGns[absPos]].append(rltvPos)
		#
		aRndmBckgrndMetrc = rtrnRndmBckgrndMetrcPrGn(aBckgrndMetrcMskd, \
		aPosRowsRndmGns,aPosClmsRndmGns)
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = rtrnRndmBckgrndMirnMetrcPrGn \
			(aBckgrndMirnMetrcMskd,aPosRowsRndmGns,aPosClmsRndmGns, \
			numMrns)
			aRndmBckgrndMirnMtrc=ma.masked_invalid(aRndmBckgrndMirnMtrc)
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Run randomization on rows
########################################################
def rnRowRndmztn(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosClmns=None,aPosProbRows=None,lenTrgtClms=False, \
	lenTrgtRows=False,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed.  aPosProbRows is an array with the 
	probabilities of sampling the positions in rows. lenTrgtRows is the 
	size of the samples to take from the rows (normarlly the same as the 
	target).aPosClmns is an array of positions in columns to be sampled. 
	aPosClmnsToGns and aPosRowsToGns are always None. Optionally, 
	aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) with an 
	array of metric values for each column and row to be summarized 
	using the value statstc. 
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of rows (of size lenTrgtRows) and values in aPosClmns
	for numRndmztns randomizations. aBcrkngdMirnStatstc is an 
	with the miRNA metrics summarized, is None if aBckgrndMirnMetrcMskd 
	is None.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbRows = len(aPosProbRows)
	if aPosClmns is None:
		aPosClmns = xrange(aBckgrndMetrcMskd.shape[1])
		try:
			raise exceptions.CelleryWarningObjct \
			('Rows are going to be randomized but no positions were',
			'provided for columns, all of them are being included.')
		except exceptions.CelleryWarningObjct as err:
			print err
			pass
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nRows==lenAPosProbRows
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		#
		aPosRows = xrange(lenAPosProbRows)
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosRowsRndm = random.choice(aPosRows,lenTrgtRows, \
		p=aPosProbRows)
		aRndmBckgrndMetrc = aBckgrndMetrcMskd[:,aPosClms]
		aRndmBckgrndMetrc = aRndmBckgrndMetrc[aPosRowsRndm,:]
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = aBckgrndMirnMetrcMskd[:,aPosClms]
			aRndmBckgrndMirnMtrc = aRndmBckgrndMirnMtrc[aPosRowsRndm,:]
			aRndmBckgrndMirnMtrc.fill_value = nan
			aRndmBckgrndMirnMtrc = ma.masked_invalid \
			(aRndmBckgrndMirnMtrc.filled())
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Run randomization on only rows and average by gene
########################################################
def rnRowRndmztnByGns(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aPosClmns=None,aPosProbRows=None,lenTrgtClms=False,lenTrgtRows=False, \
	aPosClmnsToGns=None,aPosRowsToGns=None,aBckgrndMirnMetrcMskd=None, \
	vrbse=True):
	"""
	Input: aBckgrndMetrcMskd is a 2x2 masked array with a metric of 
	interest in each cell and rows and columns in the same order as the 
	other input files. numRndmztns is the number of randomizations. 
	statstc is the statistic to calculate from the background samples. 
	seedAdd is the seed. aPosClmns is an array of positions in columns 
	to be sampled. aPosProbRows is an array with the probabilities of 
	sampling the positions in rows. lenTrgtClms is the size from the 
	column samples to take (normarlly the same as the target).
	lenTrgtRows is the size of the samples to take from the rows 
	(normally the same as the target).aPosClmnsToGns is an array of 
	position of gene to which each column in aBckgrndMetrcMskd and 
	aPosProbClmns is mapped. aPosRowsToGns is an array of position of 
	gene to which each row in aBckgrndMetrcMskd and aPosProbRows is 
	mapped. Optionally, aBckgrndMirnMetrcMskd is a masked array (3D 
	supermatrix) with an array of metric values for each column and row 
	to be summarized using the value statstc.
	Output: aBcrkngdStatstc is the statistic value for each 
	randomization of rows (of size lenTrgtRows respectively) in 
	numRndmztns. aBcrkngdMirnStatstc is an with the miRNA metrics 
	summarized, is None if aBckgrndMirnMetrcMskd is None.
	"""
	aBcrkngdStatstc = zeros(numRndmztns,dtype=float32)
	aBcrkngdStatstc.fill(nan)
	lenAPosProbRows = len(aPosProbRows)
	lenClmnsToGns = max(aPosClmnsToGns)
	lenRowsToGns = max(aPosRowsToGns)
	#vectorize the function
	rtrnRndmBckgrndMetrcPrGn =vectorize(rtrnAvrgMskdArray,excluded={0})
	if aPosClmns is None:
		aPosClmns = xrange(aBckgrndMetrcMskd.shape[1])
		try:
			raise exceptions.CelleryWarningObjct \
			('Rows are going to be randomized but no positions were',
			'provided for columns, all of them are being included.')
		except exceptions.CelleryWarningObjct as err:
			print err
			pass
	#to randomize miRNA metrics
	if aBckgrndMirnMetrcMskd is not None:
		nRows,nClmns,numMrns = aBckgrndMirnMetrcMskd.shape
		assert nRows==lenAPosProbRows
		aBcrkngdMirnStatstc = zeros((numRndmztns,numMrns),dtype=float32)
		aBcrkngdMirnStatstc.fill(nan)
		rtrnRndmBckgrndMirnMetrcPrGn = \
		vectorize(rtrnStstMskdMirnArray,excluded={0})
	else:
		aBcrkngdMirnStatstc = None
	#----------------------------
	#Run randomization	
	for smpl in xrange(numRndmztns):
		if vrbse and smpl%500 == 0:
			print '\t...Running randomization %s out of %s'%(smpl, \
			numRndmztns)
		#
		aPosClmsRndmGns = zeros(lenClmnsToGns,dytpe=float32)
		aPosClmsRndmGns.fill(nan)
		aPosRowsRndmGns = zeros(lenRowsToGns,dytpe=float32)
		aPosRowsRndmGns.fill(nan)
		#~ 
		seed('%s'%(smpl*seedAdd))#set a seed to make it 
		#replicable
		aPosRows = xrange(lenAPosProbRows)
		aPosRowsRndm = random.choice(aPosRows,lenTrgtRows, \
		p=aPosProbRows)
		for rltvPos,absPos in enumerate(aPosRowsRndm):
			aPosRowsRndmGns[aPosRowsToGns[absPos]].append(rltvPos)
		for rltvPos,absPos in enumerate(aPosClmns):
			aPosClmsRndmGns[aPosClmnsToGns[absPos]].append(rltvPos)
		#
		aRndmBckgrndMetrc = rtrnRndmBckgrndMetrcPrGn(aBckgrndMetrcMskd, \
		aPosRowsRndmGns,aPosClmsRndmGns)
		aRndmBckgrndMetrc = ma.masked_invalid(aRndmBckgrndMetrc)
		aBcrkngdStatstc[smpl] = getattr(ma,statstc)(aRndmBckgrndMetrc)
		#----------------------------
		#Run randomization miRNA metrics
		if aBckgrndMirnMetrcMskd is not None:
			aRndmBckgrndMirnMtrc = rtrnRndmBckgrndMirnMetrcPrGn \
			(aBckgrndMirnMetrcMskd,aPosRowsRndmGns,aPosClmsRndmGns, \
			numMrns)
			aRndmBckgrndMirnMtrc=ma.masked_invalid(aRndmBckgrndMirnMtrc)
			aBcrkngdMirnStatstc[smpl]=rtrnMirnStat(aRndmBckgrndMirnMtrc)
	return aBcrkngdStatstc,aBcrkngdMirnStatstc


########################################################
#~ Core method to distribute and run different randomizations
########################################################
def rndmzCore(aBckgrndMetrcMskd,mthdBckgrndRndmztn,statstcTrgt,stdTrgt, \
	outPltFl = False,aBckgrndPosProbClmnsORaPosClmns = None, \
	aBckgrndPosProbRowsORaPosRows = None, lenTrgtClms = False, \
	lenTrgtRows = False, seedAdd = False, numRndmztns = 1000, \
	maxRndmztns = 25, statstc = 'mean', fnctn = 'sf',vrbse = True, \
	aPosClmnsToGns = None, aPosRowsToGns = None, aBckgrndMirnMetrcMskd \
	= None, aTrgtMirnStatstcMskd = None, mirnDtype = 'cnt', \
	outMirnStstFl = None, aMirnNms = None):
	"""
	Input: aBckgrndMetrcMskd is masked array with background metric 
	values to be randomly sampled. mthdBckgrndRndmztn is the method to 
	make the randomization: {rnClmNRowRndmztn, rnRowRndmztn, or 
	rnClmRndmztn}. statstcTrgt is the statistic value of the target 
	whose probability is going to be calculated from the randomized 
	background using a z-score approach. stdTrgt is the standard 
	deviation of the target data. Optionally, outPltFl is a file to plot 
	the randomization and significance of target statistic. 
	aBckgrndPosProbClmnsORaPosClmns is the position of columns of 
	interest OR position probability of all columns in the background. 
	aBckgrndPosProbRowsORaPosRows is the position of rows of interest OR 
	position probability of all rows in the background. lenTrgtClms is 
	the number of columns to be sampled. lenTrgtRows is the number of 
	rows to be sampled. seedAdd is the seed to run the randomizations. 
	numRndmztns is the number of randomizations to run in the 
	background. maxRndmztns is the maximum number of iterations to 
	enforce normality. statstc is the statistic	to sample from each 
	randomization to build the normal distribution to test the 
	significance of the target statistic. fnctn is the function to 
	compare the target statistic to the resulting normal distrubtion 
	from the randomized background (sf == survival function/ greater 
	equal than). If vrbse all log messages are going to be printed.
	aPosClmnsToGns is an array of position of gene to which each column 
	in aBckgrndMetrcMskd and aPosProbClmns is mapped. aPosRowsToGns is 
	an array of position of gene to which each row in aBckgrndMetrcMskd 
	and aPosProbRows is mapped. If aPosRowsToGns and aPosClmnsToGns are 
	not None calculations are going to be run by gene. 
	aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) with an 
	array of metric values for each column and row to be summarized 
	using the value statstc. aBcrkngdMirnStatstc is an array with the 
	miRNA metrics summarized for the background. aTrgtMirnStatstcMskd is 
	an array with the miRNA metrics for the target dataset. mirnDtype is 
	the datatype of the miRNA metric: {'cnt' for counts and 'scr' for 
	scores}. outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names. fnctn
	is the statistical function to calculate the probability {'sf' is 
	greater or equal than}.
	Output: p_val is the one-side probability (following fnctn input) 
	that the target statistic belongs to the normal distributed
	statistic obtained from the randomized background, with value 
	zscore.	meanBckgrnd is the mean metric of the normal distributed 
	statistic obtained from the randomized background and stdBckgrnd is 
	its	standard error.
	"""	
	#----------------------------
	#Set required parameteres
	logRun = []
	if not seedAdd:
		seedAdd = random.random()
		try:
			raise exceptions.CelleryWarningObjct \
			('Seed for randomizations was set to',seedAdd)
		except exceptions.CelleryWarningObjct as mssge:
			if vrbse:
				print mssge
			logRun.append(mssge)
			pass
	#----------------------------
	#Run randomizations
	aBcrkngdStatstc,aBcrkngdMirnStatstc = \
	mthdBckgrndRndmztn(aBckgrndMetrcMskd,numRndmztns,statstc,seedAdd, \
	aBckgrndPosProbClmnsORaPosClmns,aBckgrndPosProbRowsORaPosRows, \
	lenTrgtClms,lenTrgtRows,aPosClmnsToGns,aPosRowsToGns, \
	aBckgrndMirnMetrcMskd,vrbse)
	#----------------------------
	#Test normality in the background randomization
	k2,pval = mstats.normaltest(aBcrkngdStatstc)
	if pval <= 0.05:
		if maxRndmztns==0:
			try:
				raise exceptions.CelleryWarningObjct \
					('Enforce background normality will stop.', \
					'The maximum number of trials was exceeded.')
			except exceptions.CelleryWarningObjct as err:
				print err
		else:
			maxRndmztns-=1
			try:
				raise exceptions.CelleryWarningObjct \
				('Re-running randomizations to enforce background', \
				'normality.')
			except exceptions.CelleryWarningObjct as mssge:
				print mssge
				pass
			if vrbse:
				print mssge
			logRun.append(mssge)
			seedAdd += numRndmztns
			p_val,zscore,meanBckgrnd,stdBckgrnd,logRun = rndmzCore \
			(aBckgrndMetrcMskd,mthdBckgrndRndmztn,statstcTrgt,stdTrgt, \
			outPltFl,aBckgrndPosProbClmnsORaPosClmns, \
			aBckgrndPosProbRowsORaPosRows,lenTrgtClms,lenTrgtRows, \
			seedAdd,numRndmztns,maxRndmztns,statstc,fnctn,vrbse, \
			aPosClmnsToGns,aPosRowsToGns,aBckgrndMirnMetrcMskd, \
			aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl,aMirnNms)
			return p_val,zscore,meanBckgrnd,stdBckgrnd,logRun
	else:
		mssge = '\n'.join([
		'\t-----------------------------------------------------------',
		'\t   The normality of the background distribution to test the',
		'\t   mean score significance had a k2: %s '%k2,
		'\t   and p-value: %s'%pval,
		'\t-----------------------------------------------------------'])
		if vrbse:
			print mssge
		logRun.append(mssge)
		#----------------------------
		#Calculate significance
		zscore = mstats.zmap(statstcTrgt,aBcrkngdStatstc)
		p_val = getattr(norm,fnctn)(zscore)#equal or greater than if 'sf'
		meanBckgrnd = np.mean(aBcrkngdStatstc)
		stdBckgrnd = np.std(aBcrkngdStatstc)
		#----------------------------
		#Make plots
		if outPltFl:
			mkPlt(statstcTrgt,aBcrkngdStatstc,outPltFl)
		#----------------------------
		#Process miRNA results
		if aBcrkngdMirnStatstc is not None:
			mssge = procMirnRndmztn(aBcrkngdMirnStatstc, \
			aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl,aMirnNms)
			if vrbse:
				print mssge
				logRun.append(mssge)
		return p_val,zscore,meanBckgrnd,stdBckgrnd,logRun


########################################################
#~ Function to be vectorized to calculate average values for matrix 
#(i.e. per gene) given a supermatrix and a set of positions
########################################################
def rtrnAvrgMskdArray(aBckgrndMetrcMskd,aRowsPosGns,aClmnsPosGns):
	"""
	Input: aBckgrndMetrcMskd is a masked array (the supermatrix) with 
	background metric values to be averaged following aRowsPosGns and 
	aClmnsPosGns. aRowsPosGns is an array with the position of rows in 
	the supermatrix	(aBckgrndMetrcMskd). aClmnsPosGns is an array with
	the positions of columns in the supermatrix.
	Output: avrgMskdVal is the average of the masked substracted 
	positions in aRowsPosGns and aClmnsPosGns of the supermatrix 
	(aBckgrndMetrcMskd).
	"""
	if len(aRowsPosGns) and len(aClmnsPosGns):
		mskdECPsPerGn = ma.mean(aBckgrndMetrcMskd[aRowsPosGns,:] \
		[:,aClmnsPosGns]).__float__()
		avrgMskdVal = float32(mskdECPsPerGn)
		return avrgMskdVal
	else:
		return nan


########################################################
#~ Calculate summary statistic for a miRNA-metric matrix.
########################################################
def rtrnMirnStat(aMirnMetrcMskd,mirnStatstc='sum'):
	"""
	Input: aMirnMetrcMskd is a masked array (3D supermatrix) with an 
	array of metric values for each column and row to be summarized 
	using the value statstc. Optionally, mirnStatstc is an statistic to 
	summarized the array of metric values.
	Output: summrzdMirnMetrcMskd is an array of the size of the input 
	array of metric values with the values summarized for all columns-
	rows.
	NOTE: the mask in aMirnMetrcMskd is a mask on column-rows of size 
	rows x columns. Masked positions are going to be excluded from the
	final summary.
	NOTE: Invalid positions are going to return nan values.
	"""
	summrzdMirnMetrcMskd = getattr(ma,mirnStatstc) \
	(ma.vstack(aMirnMetrcMskd),axis=0,dtype=float32)
	summrzdMirnMetrcMskd.fill_value = nan
	summrzdMirnMetrcMskd.filled()
	return summrzdMirnMetrcMskd


########################################################
#~ Function to be vectorized to calculate an statistic of interest for 
#RNA metrics given a supermatrix and a set of positions
########################################################
def rtrnStstMskdMirnArray(aBckgrndMirnMetrcMskd,aRowsPosGns,aClmnsPosGns, \
	lenMirnMtrc):
	"""
	Input: aBckgrndMirnMetrcMskd is a masked array (the supermatrix) with 
	stats for the metrics of shared miRNAs for each masked row-column 
	pair to be averaged following aRowsPosGns and aClmnsPosGns. 
	aRowsPosGns is an array with the position of rows in the supermatrix
	(aBckgrndMirnMetrcMskd). aClmnsPosGns is an array with the positions of 
	columns in the supermatrix. lenMirnMtrc is the length of the miRNA
	metrics array.
	Output: ststcMirnMetrcMskdVal is an array with the statistic of 
	interest caculated for aBckgrndMirnMetrcMskd following aRowsPosGns and 
	aClmnsPosGns.
	NOTE: Invalid miRNA metric values are going to return nan for all 
	positions in lenMirnMtrc.
	"""
	if len(aRowsPosGns) and len(aClmnsPosGns):
		ststcMirnMetrcMskdVal = ma.mean(aBckgrndMirnMetrcMskd[aRowsPosGns,:] \
		[:,aClmnsPosGns],axis=0,dtype=np.float32)
		ststcMirnMetrcMskdVal.fill_value = nan
		return ststcMirnMetrcMskdVal.filled()
	else:
		return array([nan for x in xrange(lenMirnMtrc)])


########################################################
#~ Calculate input length probabilities given an input model and params.
########################################################
def rtrndStrtEndCnt(lDtLens,mdlWprms,intrvlLgth,intvlJmp):
	"""
	Input: lDtLens is a list with gene/lncrna lengths the position as 
	the index in the object list.  mdlWprms is the model an parameters 
	for the distribution of lengths in lDtLens. intrvlLgth is the length 
	of the intervals (in nt) to sample, intvlJmp is the size of the 
	length interval (in nt) to increase in case of error (i.e. p == 
	inf). intrvlLgth is the length of the intervals (in nt) to sample. 
	intvlJmp is the size of the length interval to increase in case of 
	error (i.e. p == inf).
	Output: aDtPosProb is an array of genes/lncrna probabilities to 
	be sampled in the ordered of values dDtLenlPos.
	"""
	#----------------------------
	#Index lengths
	gnPos = -1
	dDtLenlPos = {}
	for gnLen in lDtLens:
		gnPos+=1
		if dDtLenlPos.has_key(gnLen):
			dDtLenlPos[gnLen].append(gnPos)
		else:
			dDtLenlPos[gnLen]=[gnPos]
	#----------------------------
	#Calculate probabilities of intervals (uniform within them)
	srtdlGnLenlPos = sorted(dDtLenlPos.keys())
	maxLength = max(srtdlGnLenlPos)
	minLength = min(srtdlGnLenlPos)
	dStrtlngthEndlngthProb = {}
	for strtLgnth in range(minLength,maxLength,intrvlLgth):
		endLgnth = strtLgnth + intrvlLgth
		if endLgnth>(maxLength):
			endLgnth = maxLength
		probd, abserr = integrate.quad(mdlWprms.pdf,strtLgnth,endLgnth)
		if probd == inf:
			intrvlLgth += intvlJmp
			try:
				raise exceptions.CelleryWarningObjct \
				('Length interval size is going to be increase by', \
				intvlJmp)
			except exceptions.CelleryWarningObjct as err:
				print err
				pass
			#increase the size of the jump
			return rtrndStrtEndCnt(lDtLens,mdlWprms,intrvlLgth,intvlJmp)
		else:
			dStrtlngthEndlngthProb[strtLgnth,endLgnth] = float32(probd)
	#----------------------------
	#Calculate probabilities for input length
	aDtPosProb = zeros(len(lDtLens),dtype=float32)#out probabilities
	sStrtlngthEndlngth = sorted(dStrtlngthEndlngthProb.keys())
	cStrt,cEnd = sStrtlngthEndlngth.pop()
	cProb = dStrtlngthEndlngthProb.pop((cStrt,cEnd))	
	cEnd+=1#set a starter
	gnLenStrt = True#set a starter
	cPosIntrvlLens = []
	while srtdlGnLenlPos:#from top to bottom
		if gnLenStrt:
			gnLen = srtdlGnLenlPos.pop()
			gnLenStrt = False
		while cStrt<=gnLen<cEnd:
			cPosIntrvlLens.extend(dDtLenlPos[gnLen])
			if srtdlGnLenlPos:
				gnLen = srtdlGnLenlPos.pop()
			else:
				gnLen = -inf#set a dummny value to pass assertion
		assert gnLen<cStrt
		if cPosIntrvlLens:
			indvldProb = divide(cProb,len(cPosIntrvlLens))
			for pos in cPosIntrvlLens:
				aDtPosProb[pos] = indvldProb
		if gnLen<0:
			break
		else:
			cStrt,cEnd = sStrtlngthEndlngth.pop()
			cProb = dStrtlngthEndlngthProb.pop((cStrt,cEnd))
			cPosIntrvlLens = []
	#----------------------------
	#test and correct for probability to sum 1	
	sumIndvldProb = np.sum(aDtPosProb,dtype=float32)
	while sumIndvldProb < float32(1):
		try:
			raise exceptions.CelleryWarningObjct \
			('probabilities were corrected to sum 1.0 from', \
			sumIndvldProb)
		except exceptions.CelleryWarningObjct as err:
			print err
			pass
		fctr = divide(1,sumIndvldProb,dtype=float32)
		aDtPosProb = multiply(aDtPosProb,fctr,dtype=float32)
		sumIndvldProb = np.sum(aDtPosProb,dtype=float32)
	return aDtPosProb


########################################################
#~ Wrapper for full randomization
########################################################
def wrprFllRndmztn(aBckgrndMetrcMskd,lBckgrndRowLens,lBckgrndClmnLens, \
	mdlWprmsRows,mdlWprmsClmns,mskRowClmnDflt=None,intrvlLgth=15, \
	intvlJmp=5,aTrgtMetrcMskd=None,statstcTrgt=False,stdTrgt=False, \
	lenTrgtRows=False,lenTrgtClms=False,outPltFl=False,numRndmztns=1000, \
	maxRndmztns=100,seedAdd=False,statstc='mean',vrbse=True, \
	aPosClmnsToGns=None,aPosRowsToGns=None,aBckgrndMirnMetrcMskd=None, \
	aTrgtMirnStatstcMskd=None,mirnDtype='cnt',outMirnStstFl=None, \
	aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. lBckgrndClmnLens is a list of 
	lengths for the rows of aBckgrndMetrcMskd. mdlWprmsRows is the 
	frozen model with parameters to sample lengths in rows 
	(==lBckgrndRowLens). mdlWprmsClmns is the frozen model with 
	parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, mskRowClmnDflt is an aditional mask for 
	aBckgrndMetrcMskd (i.e. for values that are going to be excluded 
	from calculations). intrvlLgth is an interval size to bin the 
	lengths from the samples. intvlJmp is an integer value to increase
	the size of the bins in case the porbability of one of them is inf. 
	aTrgtMetrcMskd is a masked array with metric values of interest for 
	the target. statstcTrgt is the statistic value of the 
	target whose probability is going to be calculated from the 
	randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled. lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum
	number of iterations to enforce normality. seedAdd is the seed
	to run the randomizations. statstc is the statistic to sample from 
	each randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. IfaPosRowsToGns 
	and aPosClmnsToGns are not None calculations are going to be run by 
	gene. aBckgrndMirnMetrcMskd is a masked array (3D supermatrix) with 
	an array of metric values for each column and row to be summarized 
	using the value statstc. aBcrkngdMirnStatstc is an array with the 
	miRNA metrics summarized for the background. aTrgtMirnStatstcMskd is 
	an array with the miRNA metrics for the target dataset. mirnDtype is 
	the datatype of the miRNA metric: {'cnt' for counts and 'scr' for 
	scores}. outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Both, rows and columns are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	lenClmns = len(lBckgrndClmnLens)
	assert (lenRows,lenClmns) == aBckgrndMetrcMskd.shape
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	ovrAllogRun,p_val = cmptFullRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val
	
		
########################################################
#~ Wrapper for full randomization excluding positions from background
########################################################
def wrprFllRndmztnExcldPos(aBckgrndMetrcMskd,lBckgrndRowLens, \
	lBckgrndClmnLens,mdlWprmsRows,mdlWprmsClmns,lExcldBckgrndRowPos=None, \
	lExcldBckgrndClmnPos=None,mskRowClmnDflt=None,intrvlLgth=15, \
	intvlJmp=5,aTrgtMetrcMskd=None,statstcTrgt=False,stdTrgt=False, \
	lenTrgtRows=False,lenTrgtClms=False,outPltFl=False,numRndmztns=1000, \
	maxRndmztns=100,seedAdd=False,statstc='mean',vrbse=True, \
	aPosClmnsToGns=None,aPosRowsToGns=None,aBckgrndMirnMetrcMskd=None, \
	aTrgtMirnStatstcMskd=None,mirnDtype='cnt',outMirnStstFl=None, \
	aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. lBckgrndClmnLens is a list of 
	lengths for the columns of aBckgrndMetrcMskd. mdlWprmsRows is the 
	frozen model with parameters to sample lengths in rows 
	(==lBckgrndRowLens). mdlWprmsClmns is the frozen model with 
	parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, lExcldBckgrndRowPos is the list of row positions to 
	exclude from the background. lExcldBckgrndClmnPos is the list of 
	column positions to exclude from the background. mskRowClmnDflt is 
	an aditional mask for aBckgrndMetrcMskd (i.e. for values that are 
	going to be excluded from calculations). intrvlLgth is an interval 
	size to bin the lengths from the samples. intvlJmp is an integer 
	value to increase the size of the bins in case the porbability of 
	one of them is inf. aTrgtMetrcMskd is a masked array with metric 
	values of interest for the target. statstcTrgt is the statistic 
	value of the target whose probability is going to be calculated from 
	the randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled. lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum
	number of iterations to enforce normality. seedAdd is the seed
	to run the randomizations. statstc is the statistic to sample from 
	each randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. If 
	aPosRowsToGns and aPosClmnsToGns are not None calculations are going 
	to be run by gene. aBckgrndMirnMetrcMskd is a masked array (3D 
	supermatrix) with an array of metric values for each column and row 
	to be summarized using the value statstc. aBcrkngdMirnStatstc is an 
	array with the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Both, rows and columns are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	lenClmns = len(lBckgrndClmnLens)
	assert (lenRows,lenClmns) == aBckgrndMetrcMskd.shape
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Exclude positions
	if lExcldBckgrndRowPos is not None:
		sExcldBckgrndRowPos = set(lExcldBckgrndRowPos)
		aBckgrndRowPosSlctd = []
		for pos in xrange(lenRows):
			if pos in sExcldBckgrndRowPos:
				sExcldBckgrndRowPos.remove(pos)
			else:
				aBckgrndRowPosSlctd.append(pos)
		aBckgrndRowPosSlctd = array(aBckgrndRowPosSlctd)
		lBckgrndRowLens = array(lBckgrndRowLens) \
		[aBckgrndRowPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lExcldBckgrndClmnPos is not None:
		sExcldBckgrndClmnPos = set(lExcldBckgrndClmnPos)
		aBckgrndClmnPosSlctd = []
		for pos in xrange(lenClmns):
			if pos in sExcldBckgrndClmnPos:
				sExcldBckgrndClmnPos.remove(pos)
			else:
				aBckgrndClmnPosSlctd.append(pos)
		aBckgrndClmnPosSlctd = array(aBckgrndClmnPosSlctd)
		lBckgrndClmnLens = array(lBckgrndClmnLens) \
		[aBckgrndClmnPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	ovrAllogRun,p_val = cmptFullRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val
	
		
########################################################
#~ Wrapper for full randomization including positions from background
########################################################
def wrprFllRndmztnIncldPos(aBckgrndMetrcMskd,lBckgrndRowLens, \
	lBckgrndClmnLens,mdlWprmsRows,mdlWprmsClmns,lIncldBckgrndRowPos=None, \
	lIncldBckgrndClmnPos=None,mskRowClmnDflt=None,intrvlLgth=15, \
	intvlJmp=5,aTrgtMetrcMskd=None,statstcTrgt=False,stdTrgt=False, \
	lenTrgtRows=False,lenTrgtClms=False,outPltFl=False,numRndmztns=1000, \
	maxRndmztns=100,seedAdd=False,statstc='mean',vrbse=True, \
	aPosClmnsToGns=None,aPosRowsToGns=None,aBckgrndMirnMetrcMskd=None, \
	aTrgtMirnStatstcMskd=None,mirnDtype='cnt',outMirnStstFl=None, \
	aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. lBckgrndClmnLens is a list of 
	lengths for the columns of aBckgrndMetrcMskd. mdlWprmsRows is the 
	frozen model with parameters to sample lengths in rows 
	(==lBckgrndRowLens). mdlWprmsClmns is the frozen model with 
	parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, lIncldBckgrndRowPos is the list of row positions to 
	include from the background. lIncldBckgrndClmnPos is the list of 
	column positions to include from the background. mskRowClmnDflt is 
	an aditional mask for aBckgrndMetrcMskd (i.e. for values that are 
	going to be included from calculations). intrvlLgth is an interval 
	size to bin the lengths from the samples. intvlJmp is an integer 
	value to increase the size of the bins in case the porbability of 
	one of them is inf. aTrgtMetrcMskd is a masked array with metric 
	values of interest for the target. statstcTrgt is the statistic 
	value of the target whose probability is going to be calculated from 
	the randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled. lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum
	number of iterations to enforce normality. seedAdd is the seed
	to run the randomizations. statstc is the statistic to sample from 
	each randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. If 
	aPosRowsToGns and aPosClmnsToGns are not None calculations are going 
	to be run by gene. aBckgrndMirnMetrcMskd is a masked array (3D 
	supermatrix) with an array of metric values for each column and row 
	to be summarized using the value statstc. aBcrkngdMirnStatstc is an 
	array with the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Both, rows and columns are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	lenClmns = len(lBckgrndClmnLens)
	assert (lenRows,lenClmns) == aBckgrndMetrcMskd.shape
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Include positions
	if lIncldBckgrndRowPos is not None:
		aBckgrndRowPosSlctd = array(lIncldBckgrndRowPos)
		lBckgrndRowLens = array(lBckgrndRowLens) \
		[aBckgrndRowPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lIncldBckgrndClmnPos is not None:
		aBckgrndClmnPosSlctd = array(lIncldBckgrndClmnPos)
		lBckgrndClmnLens = array(lBckgrndClmnLens) \
		[aBckgrndClmnPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	ovrAllogRun,p_val = cmptFullRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val
	
		
########################################################
#~ Wrapper for randomization on columns
########################################################
def wrprClmRndmztn(aBckgrndMetrcMskd,lBckgrndClmnLens,mdlWprmsClmns, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtClms=False, \
	outPltFl=False,numRndmztns=1000,maxRndmztns=100,seedAdd=False, \
	statstc='mean',vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None,mirnDtype='cnt', \
	outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndClmnLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsClmns is the frozen model 
	with parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, mskRowClmnDflt is an aditional maskfor aBckgrndMetrcMskd 
	(i.e. for values that are going to be excluded from calculations). 
	intrvlLgth is an interval size to bin the lengths from the samples. 
	intvlJmp is an integer value to increase the size of the bins in 
	case the porbability of one of them is inf. aTrgtMetrcMskd is a 
	masked array with metric values of interest for the target. 
	statstcTrgt is the statistic value of the target whose probability 
	is going to be calculated from the randomized background using a 
	z-score approach. stdTrgt is the standard deviation of the target 
	data. lenTrgtClms is the number of columns to be sampled. outPltFl 
	is a file to plot the randomization and significance of target 
	statistic. numRndmztns is the number of randomizations to run in the 
	background. maxRndmztns is the maximum number of iterations to 
	enforce normality. seedAdd is the seed to run the randomizations. 
	statstc is the statistic to sample from each randomization to build 
	the normal distribution to test the significance of the target 
	statistic. If vrbse all log messages are going to be printed. 
	aPosClmnsToGns is an array of position of gene to which each column 
	in aBckgrndMetrcMskd and aPosProbClmns is mapped. aPosRowsToGns is 
	an array of position of gene to which each row in aBckgrndMetrcMskd 
	and aPosProbRows is mapped. If aPosRowsToGns and aPosClmnsToGns are 
	not None calculations are going to be run by gene.
	aBckgrndMirnMetrcMskd is an array with the miRNA metrics summarized 
	for the background. aTrgtMirnStatstcMskd is an array with the miRNA 
	metrics for the target dataset. mirnDtype is the datatype of the 
	miRNA metric: {'cnt' for counts and 'scr' for scores}. outMirnStstFl 
	is the file to save the results of the miRNA probability 
	calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background. 
	NOTE: Only columns are going to be radomized.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenClmns = len(lBckgrndClmnLens)
	assert lenClmns == aBckgrndMetrcMskd.shape[1]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	lenTrgtRows = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	aBckgrndPosRows = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[0])])
	ovrAllogRun,p_val = cmptClmRndmztn(aBckgrndMetrcMskd, \
	aBckgrndPosRows,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


########################################################
#~ Wrapper for randomization on columns excluding positions from 
# background
########################################################
def wrprClmRndmztnExcldPos(aBckgrndMetrcMskd,lBckgrndClmnLens, \
	mdlWprmsClmns,lExcldBckgrndRowPos=None,lExcldBckgrndClmnPos=None, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtClms=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None,mirnDtype='cnt', \
	outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndClmnLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsClmns is the frozen model 
	with parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, lExcldBckgrndRowPos is the list of row positions to 
	exclude from the background. lExcldBckgrndClmnPos is the list of 
	column positions to exclude from the background. mskRowClmnDflt is 
	an aditional maskfor aBckgrndMetrcMskd (i.e. for values that are 
	going to be excluded from calculations). intrvlLgth is an interval 
	size to bin the lengths from the samples. intvlJmp is an integer 
	value to increase the size of the bins in case the porbability of 
	one of them is inf. aTrgtMetrcMskd is a masked array with metric 
	values of interest for the target. statstcTrgt is the statistic 
	value of the target whose probability is going to be calculated from 
	the randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtClms is the number of 
	columns to be sampled. outPltFl is a file to plot the randomization 
	and significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum 
	number of iterations to enforce normality. seedAdd is the seed to 
	run the randomizations. statstc is the statistic to sample from each 
	randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. IfaPosRowsToGns 
	and aPosClmnsToGns are not None calculations are going to be run by 
	gene. aBckgrndMirnMetrcMskd is an array with the miRNA metrics 
	summarized for the background. aTrgtMirnStatstcMskd is an array with 
	the miRNA metrics for the target dataset. mirnDtype is the datatype 
	of the miRNA metric: {'cnt' for counts and 'scr' for scores}. 
	outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only columns are going to be radomized.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenClmns = len(lBckgrndClmnLens)
	assert lenClmns == aBckgrndMetrcMskd.shape[1]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	lenTrgtRows = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
		#----------------------------
	#Exclude positions
	if lExcldBckgrndRowPos is not None:
		sExcldBckgrndRowPos = set(lExcldBckgrndRowPos)
		aBckgrndRowPosSlctd = []
		for pos in xrange(lenRows):
			if pos in sExcldBckgrndRowPos:
				sExcldBckgrndRowPos.remove(pos)
			else:
				aBckgrndRowPosSlctd.append(pos)
		aBckgrndRowPosSlctd = array(aBckgrndRowPosSlctd)
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lExcldBckgrndClmnPos is not None:
		sExcldBckgrndClmnPos = set(lExcldBckgrndClmnPos)
		aBckgrndClmnPosSlctd = []
		for pos in xrange(lenClmns):
			if pos in sExcldBckgrndClmnPos:
				sExcldBckgrndClmnPos.remove(pos)
			else:
				aBckgrndClmnPosSlctd.append(pos)
		aBckgrndClmnPosSlctd = array(aBckgrndClmnPosSlctd)
		lBckgrndClmnLens = array(lBckgrndClmnLens) \
		[aBckgrndClmnPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	aBckgrndPosRows = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[0])])
	ovrAllogRun,p_val = cmptClmRndmztn(aBckgrndMetrcMskd, \
	aBckgrndPosRows,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


########################################################
#~ Wrapper for randomization on columns including positions from 
# background
########################################################
def wrprClmRndmztnIncldPos(aBckgrndMetrcMskd,lBckgrndClmnLens, \
	mdlWprmsClmns,lIncldBckgrndRowPos=None,lIncldBckgrndClmnPos=None, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtClms=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None, \
	mirnDtype='cnt',outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndClmnLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsClmns is the frozen model 
	with parameters to sample lengths in columns (==lBckgrndClmnLens). 
	Optionally, lIncldBckgrndRowPos is the list of row positions to 
	include from the background. lIncldBckgrndClmnPos is the list of 
	column positions to include from the background. mskRowClmnDflt is 
	an aditional maskfor aBckgrndMetrcMskd (i.e. for values that are 
	going to be included from calculations). intrvlLgth is an interval 
	size to bin the lengths from the samples. intvlJmp is an integer 
	value to increase the size of the bins in case the porbability of 
	one of them is inf. aTrgtMetrcMskd is a masked array with metric 
	values of interest for the target. statstcTrgt is the statistic 
	value of the target whose probability is going to be calculated from 
	the randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtClms is the number of 
	columns to be sampled. outPltFl is a file to plot the randomization 
	and significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum 
	number of iterations to enforce normality. seedAdd is the seed to 
	run the randomizations. statstc is the statistic to sample from each 
	randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. IfaPosRowsToGns 
	and aPosClmnsToGns are not None calculations are going to be run by 
	gene. aBckgrndMirnMetrcMskd is an array with the miRNA metrics 
	summarized for the background. aTrgtMirnStatstcMskd is an array with 
	the miRNA metrics for the target dataset. mirnDtype is the datatype 
	of the miRNA metric: {'cnt' for counts and 'scr' for scores}. 
	outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only columns are going to be radomized.
	NOTE: lBckgrndClmnLens must have the same length as columns in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenClmns = len(lBckgrndClmnLens)
	assert lenClmns == aBckgrndMetrcMskd.shape[1]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	if aPosClmnsToGns is not None and aPosRowsToGns is not None:
		assert aBckgrndMetrcMskd.shape == (aPosRowsToGns,aPosClmnsToGns)
	#----------------------------
	#Define variables
	lenTrgtRows = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
		#----------------------------
	#Include positions
	if lIncldBckgrndRowPos is not None:
		aBckgrndRowPosSlctd = array(lIncldBckgrndRowPos)
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lIncldBckgrndClmnPos is not None:
		aBckgrndClmnPosSlctd = array(lIncldBckgrndClmnPos)
		lBckgrndClmnLens = array(lBckgrndClmnLens) \
		[aBckgrndClmnPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndClmnPosProb = rtrndStrtEndCnt(lBckgrndClmnLens, \
	mdlWprmsClmns,intrvlLgth,intvlJmp)
	aBckgrndPosRows = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[0])])
	ovrAllogRun,p_val = cmptClmRndmztn(aBckgrndMetrcMskd, \
	aBckgrndPosRows,aBckgrndClmnPosProb,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


########################################################
#~ Wrapper for randomization on rows
########################################################
def wrprRowRndmztn(aBckgrndMetrcMskd,lBckgrndRowLens,mdlWprmsRows, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtRows=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None, \
	mirnDtype='cnt',outMirnStstFl=None,aMirnNms=None):	
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsRows is the frozen model 
	with parameters to sample lengths in rows (==lBckgrndRowLens). 
	Optionally, mskRowClmnDflt is an aditional mask for 
	aBckgrndMetrcMskd (i.e. for values that are going to be excluded 
	from calculations). intrvlLgth is an interval size to bin the 
	lengths from the samples. intvlJmp is an integer value to increase
	the size of the bins in case the porbability of one of them is inf. 
	aTrgtMetrcMskd is a masked array with metric values of interest for 
	the target. statstcTrgt is the statistic value of the target whose 
	probability is going to be calculated from the randomized background 
	using a z-score approach. stdTrgt is the standard deviation of the 
	target data. lenTrgtRows is the number of rows to be sampled. 
	lenTrgtClms is the number of columns to be sampled. outPltFl is a 
	file to plot the randomization and significance of target statistic. 
	numRndmztns is the number of randomizations to run in the background. 
	maxRndmztns is the maximum number of iterations to enforce normality. 
	seedAdd is the seed to run the randomizations. statstc is the 
	statistic to sample from each randomization to build the normal 
	distribution to test the significance of the target statistic. If 
	vrbse all log messages are going to be printed. aPosClmnsToGns is an 
	array of position of gene to which each column in aBckgrndMetrcMskd 
	and aPosProbClmns is mapped. aPosRowsToGns is an array of position 
	of gene to which each row in aBckgrndMetrcMskd and aPosProbRows is 
	mapped. If aPosRowsToGns and aPosClmnsToGns arenot None calculations 
	are going to be run by gene. aBckgrndMirnMetrcMskd is an array with 
	the miRNA metrics summarized for the background. 
	aTrgtMirnStatstcMskd is an array with the miRNA metrics for the 
	target dataset. mirnDtype is the datatype of the miRNA metric: 
	{'cnt' for counts and 'scr' for scores}. outMirnStstFl is the file 
	to save the results of the miRNA probability calculation. aMirnNms 
	is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only rows are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	assert lenRows == aBckgrndMetrcMskd.shape[0]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	lenTrgtClms = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndPosClmns = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[1])])
	ovrAllogRun,p_val = cmptRowRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndPosClmns,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


########################################################
#~ Wrapper for randomization on rows excluding positions from background
########################################################
def wrprRowRndmztnExcldPos(aBckgrndMetrcMskd,lBckgrndRowLens, \
	mdlWprmsRows,lExcldBckgrndRowPos=None,lExcldBckgrndClmnPos=None, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtRows=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None, \
	mirnDtype='cnt',outMirnStstFl=None,aMirnNms=None):
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsRows is the frozen model 
	with parameters to sample lengths in rows (==lBckgrndRowLens). 
	lExcldBckgrndRowPos is the list of row positions to exclude from the 
	background. lExcldBckgrndClmnPos is the list of column positions to 
	exclude from the background. Optionally, mskRowClmnDflt is an 
	aditional mask for aBckgrndMetrcMskd (i.e. for values that are going 
	to be excluded from calculations). intrvlLgth is an interval size to 
	bin the lengths from the samples. intvlJmp is an integer value to 
	increase the size of the bins in case the porbability of one of them 
	is inf. aTrgtMetrcMskd is a masked array with metric values of 
	interest for the target. statstcTrgt is the statistic value of the 
	target whose probability is going to be calculated from the 
	randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled.  lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum 
	number of iterations to enforce normality. seedAdd is the seed to 
	run the randomizations. statstc is the statistic to sample from each 
	randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. If 
	aPosRowsToGns and aPosClmnsToGns are not None calculations are going 
	to be run by gene. aBckgrndMirnMetrcMskd is an array with the miRNA 
	metrics summarized for the background. aTrgtMirnStatstcMskd is an 
	array with the miRNA metrics for the target dataset. mirnDtype is 
	the datatype of the miRNA metric: {'cnt' for counts and 'scr' for 
	scores}. outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only rows are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	assert lenRows == aBckgrndMetrcMskd.shape[0]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	#----------------------------
	#Define variables
	lenTrgtClms = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Exclude positions
	if lExcldBckgrndRowPos is not None:
		sExcldBckgrndRowPos = set(lExcldBckgrndRowPos)
		aBckgrndRowPosSlctd = []
		for pos in xrange(lenRows):
			if pos in sExcldBckgrndRowPos:
				sExcldBckgrndRowPos.remove(pos)
			else:
				aBckgrndRowPosSlctd.append(pos)
		aBckgrndRowPosSlctd = array(aBckgrndRowPosSlctd)
		lBckgrndRowLens = array(lBckgrndRowLens) \
		[aBckgrndRowPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lExcldBckgrndClmnPos is not None:
		sExcldBckgrndClmnPos = set(lExcldBckgrndClmnPos)
		aBckgrndClmnPosSlctd = []
		for pos in xrange(lenClmns):
			if pos in sExcldBckgrndClmnPos:
				sExcldBckgrndClmnPos.remove(pos)
			else:
				aBckgrndClmnPosSlctd.append(pos)
		aBckgrndClmnPosSlctd = array(aBckgrndClmnPosSlctd)
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndPosClmns = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[1])])
	ovrAllogRun,p_val = cmptRowRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndPosClmns,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


########################################################
#~ Wrapper for randomization on rows including positions from background
########################################################
def wrprRowRndmztnIncldPos(aBckgrndMetrcMskd,lBckgrndRowLens, \
	mdlWprmsRows,lIncldBckgrndRowPos=None,lIncldBckgrndClmnPos=None, \
	mskRowClmnDflt=None,intrvlLgth=15,intvlJmp=5,aTrgtMetrcMskd=None, \
	statstcTrgt=False,stdTrgt=False,lenTrgtRows=False,outPltFl=False, \
	numRndmztns=1000,maxRndmztns=100,seedAdd=False,statstc='mean', \
	vrbse=True,aPosClmnsToGns=None,aPosRowsToGns=None, \
	aBckgrndMirnMetrcMskd=None,aTrgtMirnStatstcMskd=None, \
	mirnDtype='cnt',outMirnStstFl=None,aMirnNms=None):	
	"""
	Input: aBckgrndMetrcMskd is a masked array with background metric 
	values to be randomly sampled. lBckgrndRowLens is a list of lengths 
	for the rows of aBckgrndMetrcMskd. mdlWprmsRows is the frozen model 
	with parameters to sample lengths in rows (==lBckgrndRowLens). 
	lIncldBckgrndRowPos is the list of row positions to include from the 
	background. lIncldBckgrndClmnPos is the list of column positions to 
	include from the background. Optionally, mskRowClmnDflt is an 
	aditional mask for aBckgrndMetrcMskd (i.e. for values that are going 
	to be included from calculations). intrvlLgth is an interval size to 
	bin the lengths from the samples. intvlJmp is an integer value to 
	increase the size of the bins in case the porbability of one of them 
	is inf. aTrgtMetrcMskd is a masked array with metric values of 
	interest for the target. statstcTrgt is the statistic value of the 
	target whose probability is going to be calculated from the 
	randomized background using a z-score approach. stdTrgt is the 
	standard deviation of the target data. lenTrgtRows is the number of 
	rows to be sampled.  lenTrgtClms is the number of columns to be 
	sampled. outPltFl is a file to plot the randomization and 
	significance of target statistic. numRndmztns is the number of 
	randomizations to run in the background. maxRndmztns is the maximum 
	number of iterations to enforce normality. seedAdd is the seed to 
	run the randomizations. statstc is the statistic to sample from each 
	randomization to build the normal distribution to test the 
	significance of the target statistic. If vrbse all log messages are 
	going to be printed. aPosClmnsToGns is an array of position of gene 
	to which each column in aBckgrndMetrcMskd and aPosProbClmns is 
	mapped. aPosRowsToGns is an array of position of gene to which each 
	row in aBckgrndMetrcMskd and aPosProbRows is mapped. IfaPosRowsToGns 
	and aPosClmnsToGns are not None calculations are going to be run by 
	gene. aBckgrndMirnMetrcMskd is an array with the miRNA metrics 
	summarized for the background. aTrgtMirnStatstcMskd is an array 
	with the miRNA metrics for the target dataset. mirnDtype is the 
	datatype of the miRNA metric: {'cnt' for counts and 'scr' for 
	scores}. outMirnStstFl is the file to save the results of the miRNA 
	probability calculation. aMirnNms is an array of miRNA names.
	Output: ovrAllogRun is the log message of the randomization runs.
	p_val is the one-side probability (following fnctn input) that the 
	target statistic belongs to the normal distributed statistic 
	obtained from the randomized background.
	NOTE: Only rows are going to be radomized.
	NOTE: lBckgrndRowLens must have the same length as rows in 
	aBckgrndMetrcMskd.
	NOTE: mskRowClmnDflt mus have the same size as aBckgrndMetrcMskd.
	"""	
	#----------------------------
	#Test for inputs
	lenRows = len(lBckgrndRowLens)
	assert lenRows == aBckgrndMetrcMskd.shape[0]
	if mskRowClmnDflt is not None:
		assert aBckgrndMetrcMskd.shape == mskRowClmnDflt.shape
	if aPosClmnsToGns is not None and aPosRowsToGns is not None:
		assert aBckgrndMetrcMskd.shape == (aPosRowsToGns,aPosClmnsToGns)
	#----------------------------
	#Define variables
	lenTrgtClms = False
	if mskRowClmnDflt is not None:
		mskRowClmn = ma.mask_or(mskRowClmnDflt,aBckgrndMetrcMskd.mask)
		aBckgrndMetrcMskd.mask = mskRowClmn
	#----------------------------
	#Include positions
	if lIncldBckgrndRowPos is not None:
		aBckgrndRowPosSlctd = array(lIncldBckgrndRowPos)
		lBckgrndRowLens = array(lBckgrndRowLens) \
		[aBckgrndRowPosSlctd].tolist()
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[aBckgrndRowPosSlctd,:]
	if lIncldBckgrndClmnPos is not None:
		aBckgrndClmnPosSlctd = array(lIncldBckgrndClmnPos)
		aBckgrndMetrcMskd = aBckgrndMetrcMskd[:,aBckgrndClmnPosSlctd]
	#----------------------------
	#Calculate probability for column and row positions
	aBckgrndRowsPosProb = rtrndStrtEndCnt(lBckgrndRowLens,mdlWprmsRows, \
	intrvlLgth,intvlJmp)
	aBckgrndPosClmns = array([pos for pos in xrange(aBckgrndMetrcMskd. \
	shape[1])])
	ovrAllogRun,p_val = cmptRowRndmztn(aBckgrndMetrcMskd, \
	aBckgrndRowsPosProb,aBckgrndPosClmns,aTrgtMetrcMskd,statstcTrgt, \
	stdTrgt,lenTrgtRows,lenTrgtClms,outPltFl,numRndmztns,maxRndmztns, \
	seedAdd,statstc,vrbse,aPosRowsToGns,aPosClmnsToGns, \
	aBckgrndMirnMetrcMskd,aTrgtMirnStatstcMskd,mirnDtype,outMirnStstFl, \
	aMirnNms)
	return ovrAllogRun,p_val


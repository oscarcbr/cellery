#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_classes.py part of cellery (ceRNAs linking inference)
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
Test the building classes-related methods
"""

########################################################
#~ Import libraries.
########################################################
from cellery import __classes__
from numpy import array,float32

import os


#----------------------------
# test gene
name = 'ENST00000536489'
strnd = '-1'
cmmNm = 'RPH3AL'
pos = 0
chr = '17'
aIntrvls = array([(63435,63714),(65444,65593),(69413,69527),(96901, \
97076),(169210,169340),(171062,171206),(177257,177370),(202502,202633)])
lenIntrvl = __object__.clcLen(aIntrvls)
srtdMirnaNms = ['let-7a-5p','miR-181b-5p','miR-19b-2-5p']#sorted miRNAs
aMirnaCnts = array([1,2,0],dtype=float32)#counts for TargetScan
aMirnaRNAhCnts = array([1,2,1],dtype=float32)#counts for RNA hybrid
aMirnaRNAhEngy = array([-21.3,-20.3,-20.0],dtype=float32) #energy for RNA hybrid
aMirnaMirndCnts = array([3,2,0],dtype=float32)#counts for miRanda counts
aMirnaMirndScrs = array([140.0,146.1,142.0],dtype=float32)#scores for miRanda 
aMirnaMirndEngy = array([-21.3,-18.64,-15.64],dtype=float32)#energy for miRanda 
aMirnaSVMicroCnts = array([1,2,3],dtype=float32)#counts for SVMicro
aMirnaSVMicroScrs = array([240.0,136.1,122.5],dtype=float32)#scores for miRanda 
aMirnaTrgtMnrCnts = array([1,0,3],dtype=float32)#counts for TargetMiner
aMirnaPITACnts = array([1,1,3],dtype=float32)#counts for miRanda 
aMirnaPITAScrs = array([40.0,246.1,42.0],dtype=float32)#scores for miRanda 
aMirnaMirMapScrs = array([32.0,44.1,21.0],dtype=float32)#scores for mirMap
aMirnaMirMapCnts = array([2,2,3],dtype=float32)#counts for mirMap
aMirnaWalk = array([0,1,0],dtype=bool)#mirWalk confirmation
#
gene = __classes__.gene(name,srtdMirnaNms)
gene.strnd = strnd
gene.cmmNm = cmmNm
gene.pos = pos
gene.chr = chr
gene.aIntrvls = aIntrvls
gene.len = lenIntrvl
gene.aMirnaCnts = aMirnaCnts
gene.aMirnaRNAhCnts = aMirnaRNAhCnts
gene.aMirnaRNAhEngy = aMirnaRNAhEngy
gene.aMirnaMirndCnts = aMirnaMirndCnts
gene.aMirnaMirndScrs = aMirnaMirndScrs
gene.aMirnaMirndEngy = aMirnaMirndEngy
gene.aMirnaSVMicroCnts = aMirnaSVMicroCnts
gene.aMirnaSVMicroScrs = aMirnaSVMicroScrs
gene.aMirnaTrgtMnrCnts = aMirnaTrgtMnrCnts
gene.aMirnaPITACnts = aMirnaPITACnts
gene.aMirnaPITAScrs = aMirnaPITAScrs
gene.aMirnaMirMapScrs = aMirnaMirMapScrs
gene.aMirnaMirMapCnts = aMirnaMirMapCnts
gene.aMirnaWalk = aMirnaWalk


#----------------------------
# test node
w = float32(1.255)
node = __classes__.node('ENST00000000001')

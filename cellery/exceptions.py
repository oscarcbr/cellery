#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  exceptions.py part of cellery (ceRNAs linking inference)
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
Celery exceptions and packages
"""

class CelleryStandardError(StandardError):
	"""Cellery standard error.
	Cellery should use this error (or subclasses of it), making it easy 
	to distinguish all standard error messages.
	"""
	pass

class CelleryExceptionObjct(Exception):
	"""Cellery exception for the object class."""
	def __init__(self, message):
		self.message = message
	def __str__(self):
		return self.message

class CelleryWarningObjct(Exception):
	"""Cellery warning for the object class."""
	def __init__(self, message, value):
		self.value = value
		self.message = message
	def __str__(self):
		return 'CelleryWarningObjct: %s %s'%(self.message,self.value)

# -*- coding: utf-8 -*-
"""
General tools for recon3d
"""

import sys,os
import numpy as np

###===========================================================================
def pixel_offsets(p):
   """
   Calculates the pixel offsets for each scene (positive southward and eastward)

   Parameters
   ----------

   p :   Parser class object

   """
   for i in np.arange(len(p.scenes)):
      latoff = (p.gen.cornlat - p.scenes[i].cornlat)/p.scenes[i].spacelat
      p.scenes[i].offlat = -1*_round(latoff)  
      
      lonoff = (p.gen.cornlon - p.scenes[i].cornlon)/p.scenes[i].spacelon
      p.scenes[i].offlon = -1*_round(lonoff) 

      p.scenes[i].eline = p.scenes[i].rows + p.scenes[i].offlat  
      p.scenes[i].sline = -p.scenes[i].offlat 

   column_buffer(p)
###----------------------------------------------
def column_buffer(p):
   """
   Calculates column buffers

   Parameters
   ---------

   p :  Parser class object 

   """
   for i in np.arange(len(p.scenes)):
      if p.scenes[i].offlon <= 0:
         p.scenes[i].bcol = 0
         p.scenes[i].scol = -p.scenes[i].offlon   # scols is a 0-base index
      else:
         p.scenes[i].bcol = p.scenes[i].offlon 
         p.scenes[i].scol = 0

###===========================================================================
def _round(a):
   """
   Rounding routine to always round away from 0 if |a - int(a)| >= 0.5

   Parameters
   ---------
   
   a :  float
   """
   return np.int32(a + np.sign(a)*0.5)

###===========================================================================
def _usage_warning(mode=None):
   """ 
   Recon3D:  3D reconstruction from InSAR data

   Usage: 
 
   Generate command_file template: 
      recon3d_driver.py -g [number_of_scenes] > command_file
   
   @@@
   Run in serial:  
      recon3d_driver.py command_file [1]

   @@@
   Run in parallel:
      mpiexec [-n num_procs] recon3d_driver.py command_file -p

   """
   doc = _usage_warning.__doc__
   if mode == 'parallel':
      doc = doc.replace('recon3d_driver.py','recon3d_parallel.py')
      doc = doc.split('@@@')
      del doc[1]
      doc[1] = doc[1].replace('-p','')
      doc = ''.join(doc)
   elif mode == 'serial':
      doc = doc.replace('recon3d_driver.py','recon3d_serial.py')
      doc = ''.join(doc.split('@@@')[:2])
   else:
      doc = ''.join(doc.split('@@@')[:3])
   print(doc)
   sys.exit()

###===========================================================================
def _make_command_file(n):
   from _reconparser import CommandFile
   cmd = CommandFile()
   basic = cmd.basic.split('\n')
   for i in np.arange(len(basic)):
      if basic[i] != '': basic[i] = basic[i].lstrip()
   basic = '\n'.join(basic)

   scene = cmd.scene.split('\n')
   for i in np.arange(len(scene)):
      if scene[i] != '': scene[i] = scene[i].lstrip()
   scene = '\n'.join(scene)

   closer = cmd.closer.split('\n')
   for i in np.arange(len(closer)):
      if closer[i] != '': closer[i] = closer[i].lstrip()
   closer = '\n'.join(closer)

   print(basic)
   print(cmd.options)
   for i in np.arange(1,n+1):
      print(scene.replace('<<scene_num>>',str(i)))
   print(closer)
   sys.exit()

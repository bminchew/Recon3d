# -*- coding: utf-8 -*-
"""

Part of recon3d.py

Utility module to read and parse the command file 


"""

import sys,os
import datetime
import numpy as np
import _recontools

__all__ = ['Parser','Scene','General']

###===========================================================


###===========================================================


class Parser():

   """
   Class to read and parse command file

   Usage:  a = parser.Parser(command_file_name)
   """

   def __init__(self,cmd):
      fid = open(cmd,'r')
      freads = fid.readlines()
      fid.close()

      self.gen = General()
      self.scenes = []
      master = False

      for frow in freads:
         r = frow.split('::')
         if len(r) > 1:
            r[0],r[1] = r[0].strip(),r[1].strip()
            if 'Number of Scenes (-)' == r[0]:
               self.gen.numscenes = np.int32(r[1])
            elif 'Master Scene (-)' == r[0]:
               self.gen.masterscn = np.int32(r[1])
            elif 'Radar Wavelength (in preferred output distance units)' == r[0]:
               self.gen.lamb = np.float32(r[1])
               self.gen.r2dsp = self.gen.lamb / (4.*np.pi)
            elif 'Time converstion factor (applied to all scenes equally)' == r[0]:
               if 'n' in r[1] or 'N' in r[1]:
                  self.gen.timefac = 1.
               else:
                  self.gen.timefac = np.float32(r[1])
            elif 'Distance converstion factor (applied to all scenes equally)' == r[0]:
               if 'n' in r[1] or 'N' in r[1]:
                  self.gen.distfac = 1.
               else:
                  self.gen.distfac = np.float32(r[1])

            elif 'Prior model covariance weighting factor' == r[0]:
               self.gen.cmlambda = np.float32(r[1])
            elif 'lsqr damping factor' == r[0]:
               self.gen.damp = np.float32(r[1])

            elif 'Output folder' in r[0]:
               self.gen.outfldr = r[1]
               if self.gen.outfldr[-1] != '/': 
                  self.gen.outfldr += '/'
            elif 'Output file prefix' in r[0]:
               self.gen.outpref = r[1]
            elif 'Output null value' in r[0]:
               self.gen.outnull = np.float32(r[1])


            elif 'Output diag(G^T*W*G)^-1)' in r[0]:
               self._bool_dict('Output diag(G^T*W*G)^-1)',r[1])
            elif 'Output diag(G^T*G)^-1' in r[0]:
               self._bool_dict('Output diag(G^T*G)^-1',r[1])
            elif 'Output offdiag(G^T*G)^-1' in r[0]:
               self._bool_dict('Output offdiag(G^T*G)^-1',r[1])
            elif 'Output offdiag(G^T*W*G)^-1)' in r[0]:
               self._bool_dict('Output offdiag(G^T*W*G)^-1)',r[1])
            elif 'Output GDOP sqrt(tr((G^TG)^-1))' in r[0]:
               self._bool_dict('Output GDOP sqrt(tr((G^TG)^-1))',r[1])
            elif 'Output sqrt(tr((G^T*W*G)^-1))' in r[0]:
               self._bool_dict('Output sqrt(tr((G^T*W*G)^-1))',r[1])
            elif 'Output mean square error estimate' in r[0]:
               self._bool_dict('Output mean square error estimate',r[1])
            elif 'Output velocity magnitude' in r[0]:
               self._bool_dict('Output velocity magnitude',r[1])
            elif 'Output average correlation' in r[0]:
               self._bool_dict('Output average correlation',r[1])
            elif 'Output number of scenes' in r[0]:
               self._bool_dict('Output number of scenes',r[1])
            elif 'Output rank(G) at each pixel' in r[0]:
               self._bool_dict('Output rank(G) at each pixel',r[1])

            elif 'Scene' == r[0]:
               self.scenes.append(Scene())
               if np.int32(r[1]) == self.gen.masterscn:
                  self.gen.masterind = len(self.scenes) - 1
            elif 'Unwrapped or LOS displacement file' in r[0]:
               self.scenes[-1].unwfile = r[1]
            elif 'Samples in Unwrapped or LOS displacement (-)' in r[0]:
               self.scenes[-1].cols = np.int32(r[1])
            elif 'Type of unwrapped or LOS displacement file' in r[0]:
               if r[1] in self.scenes[-1].disopts:
                  self.scenes[-1].disphz = 'disp'
               else: 
                  self.scenes[-1].disphz = 'phase'
            elif 'Phase or displacement null value' in r[0]:
               self.scenes[-1].phznull = np.float32(r[1])
            elif 'Correlation file' in r[0]:
               self.scenes[-1].corfile = r[1]
            elif 'Correlation null value' in r[0]:
               self.scenes[-1].cornull = np.float32(r[1])
            elif 'LOS file' in r[0]:
               self.scenes[-1].losfile = r[1]
            elif 'LOS null value' in r[0]:
               self.scenes[-1].losnull = np.float32(r[1])
            elif 'Upper left corner Latitude' in r[0]:
               self.scenes[-1].cornlat = np.float64(r[1])
            elif 'Upper left corner Longitude' in r[0]:
               self.scenes[-1].cornlon = np.float64(r[1])
            elif 'Latitude Spacing' in r[0]:
               self.scenes[-1].spacelat = np.float64(r[1])
            elif 'Longitude Spacing' in r[0]:
               self.scenes[-1].spacelon = np.float64(r[1])
            elif 'Number of looks' in r[0]:
               self.scenes[-1].dlooks = 2.*np.float32(r[1])

            elif 'End of command file' == r[0]:
               break
      if len(self.scenes) != self.gen.numscenes:
         ValueError('Number of scenes does not match intended number of scenes')

      self._resolve_conflicts()
      self._master_values()
      self._get_rows()
      self._pixel_offsets()
      self._disp2unw()


   def _resolve_conflicts(self):
      """
      Reconciles any values that conflict
      """
      if np.abs(self.gen.cmlambda) > 1.e-7:  self.gen.damp = 0.


   def _master_values(self):
      """
      Get the number of columns, corner position, and grid spacing for the master scene
      """
      self.gen.cols     = self.scenes[self.gen.masterind].cols
      self.gen.cornlat  = self.scenes[self.gen.masterind].cornlat
      self.gen.cornlon  = self.scenes[self.gen.masterind].cornlon
      self.gen.spacelat = self.scenes[self.gen.masterind].spacelat
      self.gen.spacelon = self.scenes[self.gen.masterind].spacelon 
      self.gen._load_grdarray()

   def _get_rows(self):
      """
      Calculate the number of rows in each scene and find the length of the output image.  

      Calls Scene._get_rows for each scene 
      """
      rows = np.empty(len(self.scenes),dtype=np.int32)
      for i in np.arange(len(self.scenes)):
         self.scenes[i]._get_rows()
         rows[i] = self.scenes[i].rows
      rows = np.sort(rows)
      self.gen.erow = rows[-3]
      self.gen.rows = self.scenes[self.gen.masterind].rows
      if self.gen.rows  <= self.gen.erow:
         self.gen.doline = self.gen.rows
      else:
         self.gen.doline = self.gen.erow
      
   def _pixel_offsets(self):
      """
      Calculates the pixel offset for each scene (positive southward and eastward)
      """
      for i in np.arange(len(self.scenes)):
         self.scenes[i]._get_offsets(self.gen.grd)
         self.scenes[i]._get_col_buff() 

   def _disp2unw(self):
      """
      Convert displacement to phase
      """
      for i in np.arange(len(self.scenes)):
         self.scenes[i].dsp2r = 1. 
         if self.scenes[i].disphz == 'disp':
            self.scenes[i].dsp2r *= 4.*np.pi / self.gen.lamb 

   def _bool_dict(self,challenge,response):
      if 'y' in response or 'Y' in response or 't' in response or 'T' in response:
         self.gen.outputs[challenge] = True
      else:
         self.gen.outputs[challenge] = False

###===========================================================
class Scene():
   """ 
   Scene  :  Class containing information for a given scene

   Variables 
   ---------

   unwfile :   unwrapped interferogram filename
   corfile :   correlation filename
   losfile :   LOS filename (assumed to be bit-interleaved)

   cols    :   number of columns (samples) in the unwfile
   rows    :   number of rows (lines) in unwfile -- determined from filesize and cols 
   looks   :   number of looks
   
   cornlat :   upper left corner latitude (deg)
   cornlon :   upper left corner longitude (deg)
   spacelat :  latitude spacing (deg/pixel)
   spacelon :  longitude spacing (deg/pixel)

   disphz  :   LOS displacement or phase values are given in unwfile (default is phase)

   phznull :   Phase or displacement null value
   cornull :   Correlation null value
   losnull :   LOS null value
   

   Notes
   -----

   - All input files are assumed to be 32-bit floating point headerless binary grid files

   """

   def __init__(self):
      self.disphz = 'phase'
      self.disopts = ['disp','Disp','DISP'] 
      self.rows = -1

   def _get_rows(self,nbytes=4.):
      """
      Calculate number of rows from the filesize and number of columns
      """
      sz = os.path.getsize(self.unwfile)
      self.rows = np.int32( sz / (nbytes*self.cols) )

   def _get_offsets(self,mastergrd):
      """ 
      Calculates the pixel offsets (positive southward and eastward)
      """

      latoff = (mastergrd[0] - self.cornlat)/self.spacelat
      self.offlat = -1*_recontools._round(latoff)

      lonoff = (mastergrd[1] - self.cornlon)/self.spacelon
      self.offlon = -1*_recontools._round(lonoff)

      self.eline = self.rows + self.offlat - 1
      self.sline = -self.offlat

   def _get_col_buff(self):
      """
      Calculate the left hand column buffer size in pixels
      """
      if self.offlon <= 0:
         self.bcol, self.scol = 0, -self.offlon
      else:
         self.bcol, self.scol = self.offlon, 0
      

###===========================================================
class General():
   """
   Class to store general information, typically regarding the output image
   """
   def __init__(self):
      self.note = 'Created on', datetime.datetime.now()

      self.numscenes = None    #  total number of scenes
      self.masterscn = None    #  master scene

      self.masterind = None    #  index of master scene in list of scenes
      self.cols      = None    #  columns (samples) in master scene
      self.rows      = None    #  rows (lines) in master scene
      self.erow      = None    #  number of rows in the output scene
 
      self.cornlat   = None    #  corner latitude of master scene
      self.cornlon   = None    #  corner longitude of master scene
      self.spacelat  = None    #  latitude spacing of master scene
      self.spacelon  = None    #  longitude spacing of master scene
      self.grd       = np.array([self.cornlat, self.cornlon, self.spacelat, self.spacelon])

      self.lamb      = None    #  radar wavelength
      self.timefac   = None    #  time conversion factor
      self.distfac   = None    #  distance conversion factory
      self.cmlambda  = 1.e-4   #  prior model covariance weighting factor (usaully called lambda)
      self.damp      = 0.      #  parameter for standard damping in scipy.linalg.lsqr algorithm

      self.outfldr   = './'    #  output folder
      self.outpref   = 'fluye' #  output file prefix
      self.outnull   = 0       #  output null value

      self.outputs   = {}      #  dictionary to hold output options
      self.addfile   = 'fluye.address'
      self.numfile   = 'fluye.numscenes'

   def _load_grdarray(self):
      self.grd = np.array([self.cornlat, self.cornlon, self.spacelat, self.spacelon])




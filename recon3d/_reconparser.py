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

class _Options():
   """
   define the inputs and set the default values 
   """
   def __init__(self):
      self.optdefs = [('Output diag(G^T*W*G)^-1)','n'),
                  ('Output diag(G^T*G)^-1','n'),
                  ('Output offdiag(G^T*G)^-1','n'),
                  ('Output offdiag(G^T*W*G)^-1)','n'),
                  ('Output GDOP sqrt(tr((G^TG)^-1))','n'),
                  ('Output sqrt(tr((G^T*W*G)^-1))','n'),
                  ('Output mean square error estimate','n'),
                  ('Output velocity magnitude','n'),
                  ('Output average correlation','n'),
                  ('Output number of scenes','y'),
                  ('Output rank(G) at each pixel','n')]

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
      self._load_defaults_bool_dict()

      for frow in freads:
         r = frow.split('::')
         if len(r) > 1:
            r = [x.strip() for x in r]
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

            elif 'PETSc KSP rtol' == r[0] or 'PETSc KSP reltol' == r[0]:
               self.gen.petsc_rtol = np.float32(r[1])
            elif 'PETSc KSP atol' == r[0] or 'PETSc KSP abstol' == r[0]:
               self.gen.petsc_atol = np.float32(r[1])
            elif 'PETSc KSP dtol' == r[0] or 'PETSc KSP divtol' == r[0]:
               self.gen.petsc_dtol = np.float32(r[1])
            elif 'PETSc KSP max_it' == r[0] or 'PETSc KSP maxits' == r[0]:
               self.gen.petsc_max_it = np.float32(r[1])

            elif 'scipy.sparse.linalg atol' == r[0]:
               self.gen.scipy_atol = np.float32(r[1])
            elif 'scipy.sparse.linalg btol' == r[0]:
               self.gen.scipy_btol = np.float32(r[1])
            elif 'scipy.sparse.linalg iter_lim' == r[0]:
               self.gen.scipy_iter_lim = np.float32(r[1])
            elif 'scipy.sparse.linalg conlim' == r[0]:
               self.gen.scipy_conlim = np.float32(r[1])

            elif 'Output folder' in r[0]:
               self.gen.outfldr = r[1]
               if self.gen.outfldr[-1] != '/': 
                  self.gen.outfldr += '/'
            elif 'Output file prefix' in r[0]:
               self.gen.outpref = r[1]
            elif 'Output null value' in r[0]:
               self.gen.outnull = np.float32(r[1])

            elif r[0] in self.gen.outputs.keys():  # define options 
               self._bool_dict(r[0],r[1])         

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

   ####################################################

   def _load_defaults_bool_dict(self):
      opts = _Options()
      for entry in opts.optdefs:
         self._bool_dict(entry[0],entry[1]) 

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
      self.cmlambda  = 1.      #  prior model covariance weighting factor (usaully called lambda)
      self.damp      = 0.      #  parameter for standard damping in scipy.linalg.lsqr algorithm

      self.outfldr   = './'    #  output folder
      self.outpref   = 'r3dout'#  output file prefix
      self.outnull   = 0       #  output null value

      self.outputs   = {}      #  dictionary to hold output options
      self.addfile   = 'r3d.address'
      self.numfile   = 'r3d.numscenes'
      self.metfile   = '.r3d.metfile'
      self.metdelim  = ':'

      self.petsc_rtol   = None # relative tolerance for PETSc KSP solver (None sets to PETSc default)
      self.petsc_atol   = None # absolute tolerance for PETSc KSP solver
      self.petsc_dtol   = None # divergence tolerance for PETSc KSP solver
      self.petsc_max_it = None # max iterations for PETSc KSP solver

      self.scipy_atol   = 1.e-8  # atol for scipy.sparse.linalg.lsqr
      self.scipy_btol   = 1.e-8  # btol for scipy.sparse.linalg.lsqr
      self.scipy_conlim = 1.e8   # conlim for scipy.sparse.linalg.lsqr
      self.scipy_iter_lim = 2.e5 # max iterations for scipy.sparse.linalg.lsqr

   def _load_grdarray(self):
      self.grd = np.array([self.cornlat, self.cornlon, self.spacelat, self.spacelon])



###===============================================================

class CommandFile():
   def __init__(self):
      opt = _Options()
      
      self.basic = """
         Command file for recon3d


         Number of Scenes (-)                                        ::  
         Master Scene (-)                                            ::  1

         Radar Wavelength (in preferred output distance units)       ::  0.24
         Time converstion factor (applied to all scenes equally)     ::  365 
         Distance converstion factor (applied to all scenes equally) ::  none

         Output folder                       :: . 
         Output file prefix                  :: recon3d_output

         Output null value                   :: -100
         Prior model covariance weighting factor :: 1.e-4
         lsqr damping factor                 :: 0
         """

      self.options = '\n'
      spacer = 0
      for entry in opt.optdefs:
         if spacer % 4 == 0: self.options += '\n'
         self.options += entry[0]
         k = len(entry[0]) 
         m = np.max([39,k])
         for i in np.arange(k,m): self.options += ' '
         self.options += ' ::  ' + entry[1] + '\n'
         spacer += 1

      self.scene = """

         Scene :: <<scene_num>> 
            Unwrapped or LOS displacement file           :: 
            Samples in Unwrapped or LOS displacement (-) :: 
            Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
            Phase or displacement null value          :: 0

            Correlation file                          :: 
            Correlation null value                    :: 0

            LOS file (must be by-pixel interleave in order ENU) :: 
            LOS null value                            :: -100

            Upper left corner Latitude (deg)          :: 
            Upper left corner Longitude (deg)         :: 
            Latitude Spacing (deg/pixel)              :: 
            Longitude Spacing (deg/pixel)             :: 

            Number of looks                           :: 
               """

      self.closer = """
         End of command file

         This file was automatically generated using recon3d_driver.py -g
               """


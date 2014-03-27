#!/usr/bin/env python

"""
Recon3D

3D reconstruction of displacement (velocity) field using multiple 
InSAR acquisitions.  Data must have been acquired from at least
3 different line-of-sight directions.

for usage:  recon3d_driver.py -h
"""
import sys,os
import numpy as np
from recon3d import _reconparser, _recontools, _reconsolver

def main(args,mode):
   if mode == 0:
      # serial using scipy
      p = _reconparser.Parser(args[0])
      solve = _reconsolver.Solver(p)
      solve.define_solution_space(p)
      solve.BS(p)

   elif mode == 1:
      # parallel using PETSc
      try:
         import mpi4py
         from mpi4py import MPI 
      except:
         raise mpi4pyError('mpi4py must be installed and in PYTHONPATH')
      try:
         import petsc4py
         petsc4py.init([],comm=MPI.COMM_WORLD)
         from petsc4py import PETSc
      except:
         raise petsc4pyError('petsc4py must be installed and in PYTHONPATH')

      p = _reconparser.Parser(args[0])
      solve = _reconsolver.Solver(p)
      solve.PPS(p,comm)

class mpi4pyError(Exception):
   pass

class petsc4pyError(Exception):
   pass        

if __name__ == '__main__':
   # mode = 0 for serial, = 1 for parallel
   args, mode = sys.argv[1:], 0
   if len(args) < 1 or len(args) > 2:
      print(__doc__)
      sys.exit()
   elif '-h' in args or '-help' in args or '--help' in args:
      _recontools._usage_warning()
   elif '-g' in args:
      try: 
         nscenes = np.int32(args[1])
      except: 
         nscenes = 3
      _recontools._make_command_file(nscenes)
   elif '-p' in args:
      mode = 1
   main(args=args,mode=mode)



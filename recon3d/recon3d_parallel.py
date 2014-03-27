#!/usr/bin/env python

"""
Recon3D

3D reconstruction of displacement (velocity) field using multiple 
InSAR acquisitions.  Data must have been acquired from at least
3 different line-of-sight directions.

This is a parallelized implementation that requires PETSc, petsc4py, 
and mpi4py connected to a working MPI library 

for usage:  recon3d_parallel.py -h
"""

import sys,os
import subprocess
import numpy as np
from recon3d import _reconparser, _recontools, _reconsolver

def main(args):
   try:
      import mpi4py
      from mpi4py import MPI 
   except:
      raise mpi4pyError('mpi4py must be installed and in PYTHONPATH')
   try:
      import petsc4py
      comm = MPI.COMM_WORLD      
      petsc4py.init([],comm=comm)
      from petsc4py import PETSc
   except:
      raise petsc4pyError('petsc4py must be installed and in PYTHONPATH')

   p = _reconparser.Parser(args[0])
   solve = _reconsolver.Solver(p)
   comm.Barrier()
   solve.PPS(p,comm)

class mpi4pyError(Exception):
   pass
class petsc4pyError(Exception):
   pass

if __name__ == '__main__':
   args = sys.argv[1:]
   if len(args) < 1:
      print(__doc__)
      sys.exit()
   elif '-h' in args or '-help' in args or '--help' in args:
      _recontools._usage_warning(mode='parallel')
   elif '-g' in args[0]:
      try: 
         nscenes = np.int32(args[1])
      except: 
         nscenes = 1 
      _recontools._make_command_file(nscenes)
   main(args=args)


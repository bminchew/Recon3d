# -*- coding: utf-8 -*-
"""

Module to read data and load the design (G) and covariance (C_phi) matrices 

"""

import datetime
import sys,os,time
import numpy as np
import _recontools, _reconparser

__all__ = ['Solver']

###===========================================================

###===========================================================


class Solver():
   """
   Solver Class

   Contains all of the solver functions.  This is the main engine for the recon3d code.

   Usage:  Solver(parser_class)  
   """

   def __init__(self,p):
      sys.stdout.write('Instantiating Solver class\n')
      self.note = 'Created on', datetime.datetime.now()

###-------------------------------------------------------------
   def define_solution_space(self,p):
      from _reconutils import defsolspace
      sys.stdout.write('Defining the solution space\n')

      metfile, delim = p.gen.metfile, p.gen.metdelim
      self._write_metfile(p,filename=metfile,delim=delim)
      lp = len(p.scenes)

      address,rankout,nonulout,p.gen.gmat_rows,addind = defsolspace(metfile,lp,p.gen.doline,p.gen.cols,1)
      address = address[:addind]
      self._write_address(address=address,p=p)

      fid = open('.gmat_rows','w')
      fid.write(str(p.gen.gmat_rows))
      fid.close()

      if p.gen.outputs['Output number of scenes']:
         self._write_output(data=nonulout,p=p,ext='.num')
      if p.gen.outputs['Output rank(G) at each pixel']:
         self._write_output(data=rankout,p=p,ext='.rnk')

###-------------------------------------------------------------
   def BS(self,p):
      """

      Basic solver that inverts the entire scene using a single patch.
      The inversion is carried out with scipy.sparse.linalg.lsqr

      """

      sys.stdout.write("\nRunning BS:  Serial basic solver\n")
      from scipy.sparse import csr_matrix, dia_matrix, bsr_matrix, coo_matrix, identity
      from scipy.sparse.linalg import lsqr
      from _reconutils import get_gcd_vals, laplace_builder

      cols = p.gen.cols
      address = self._read_address(p)

      try:
         p.gen.gmat_rows
      except:
         fid = open('.gmat_rows','r')
         p.gen.gmat_rows = np.int32(''.join(fid.readlines()))
         fid.close()


      sys.stdout.write('building sparse matrices\n')

      metfile = p.gen.metfile
      gmatvals, gmatij, cdiv, dvec = get_gcd_vals(metfile,len(p.scenes),
               p.gen.doline,p.gen.cols,address,p.gen.gmat_rows)


      inds = np.arange(p.gen.gmat_rows)
      gmatij = gmatij.reshape(-1,3*p.gen.gmat_rows) 
      gmat = bsr_matrix((gmatvals,gmatij),dtype=np.float32)
      cdi  = coo_matrix((cdiv,(inds,inds)),dtype=np.float32).todia()
      del gmatvals, gmatij, cdiv, inds
      sys.stdout.write('gmat = [%8d,%8d]; nnz = %8d\n'%(gmat.get_shape()[0],gmat.get_shape()[1],gmat.getnnz()))
      sys.stdout.write('cdi = [%8d,%8d]; nnz = %8d\n'%(cdi.get_shape()[0],cdi.get_shape()[1],cdi.getnnz()))
      sys.stdout.write('building A matrix\n')
      omega = gmat.transpose().dot(cdi)
      del cdi
      dtilde = omega.dot(dvec)
      del dvec
      gtg = omega.dot(gmat)
      del omega, gmat


      sys.stdout.write('building fmat\n') 
#      fmat = self._build_laplace_scipy_sparse(a=address,p=p)
      fmatval, fmatij = laplace_builder(p.gen.cols,address)
      fmatval = fmatval[:np.argmin(fmatval)]
      fmatij = fmatij.reshape(2,-1)[:,:len(fmatval)]
      fmat = csr_matrix((fmatval,fmatij),dtype=np.float32)
      inds = np.arange(gtg.get_shape()[0])

      doobscov = True
      if doobscov:
         cli = coo_matrix((gtg.diagonal(),(inds,inds)),dtype=np.float64).todia()
         resfilename = 'ResidualNorm_diagGtG.txt'
      else:
         cli = identity(gtg.get_shape()[0],dtype=np.float64,format='dia')
         resfilename = 'ResidualNorm_identity.txt'

      rhv = fmat.transpose().dot(cli).dot(fmat)
      del inds,fmat,cli

      doLcurve = False
      if doLcurve:
         nvals = 30
         lcvals = np.array([0])
         lcvals = np.append(lcvals,np.logspace(-6,1,nvals))
         rnorms = np.empty(len(lcvals),dtype=np.float32)
         i = 0
         for p.gen.cmlambda in lcvals:
            sys.stdout.write('cmlambda = %6f\n'%p.gen.cmlambda)
            A = gtg + p.gen.cmlambda*rhv
         
            sys.stdout.write('inverting using scipy.sparse.linalg.lsqr\n')
            t0 = time.time()
            mtilde = lsqr(A,dtilde, damp=p.gen.damp, atol=p.gen.scipy_atol,
                     btol=p.gen.scipy_btol, conlim=p.gen.scipy_conlim,
                     iter_lim=p.gen.scipy_iter_lim, show=False, calc_var=False)  
            sys.stdout.write('lsqr took %10f seconds\n'%(time.time()-t0))
            sys.stdout.write('termination type = %2d with %6d iterations\n'%(mtilde[1],mtilde[2]))
            sys.stdout.write('norm(r) = %12f\n' % mtilde[3])
            rnorms[i] = mtilde[3]
            i += 1
         
         fid = open(resfilename,'w')
         fid.write('lambda values       norm(r = d - Ax)  \n')
         for i in np.arange(len(lcvals)):
            fid.write(str(lcvals[i]) + '     ' + str(rnorms[i]) + '\n')
         fid.close()
      
      else:
         if np.abs(p.gen.cmlambda) > 1.e-12:  # this can be sped up by not finding fmat unless necessary
            gtg = gtg.tocsr()
            gtg = gtg + p.gen.cmlambda*rhv

         sys.stdout.write('inverting using scipy.sparse.linalg.lsqr\n')
         t0 = time.time()
         mtilde = lsqr(gtg,dtilde,damp=p.gen.damp,conlim=1.e8,iter_lim=2e5,show=False)
         sys.stdout.write('lsqr took %10f seconds\n'%(time.time()-t0))
         sys.stdout.write('termination type = %2d with %6d iterations\n'%(mtilde[1],mtilde[2]))
         mtilde = mtilde[0]*p.gen.timefac * p.gen.distfac * p.gen.r2dsp

         self._write_post_model(model=mtilde,p=p)
         self._cleanup(p=p)
###-----------------------------------------------------------------------------


###-----------------------------------------------------------------------------
   def PPS(self,p,comm):
      """
      Parallel solver implementing a PETSc KSP linear routine over the entire scene as a 
      single patch (PPS = Parallel PETSc Single)
      """
      from _reconutils import get_gcd_vals, laplace_builder, pixelxpixelsolver
      try:
         import mpi4py
         from mpi4py import MPI 
      except:
         sys.stdout.write('mpi4py must be installed and in PYTHONPATH\n')
         sys.exit()
      try:
         import petsc4py
         petsc4py.init([],comm=MPI.COMM_WORLD)
         from petsc4py import PETSc
      except:
         sys.stdout.write('petsc4py must be installed and in PYTHONPATH\n')
         sys.exit()

      PETSc.Sys.Print('Running PPS:  Parallel PETSc (linear) Solver')

      size, rank = comm.Get_size(), comm.Get_rank()
      cols = p.gen.cols

      if rank == 0:
         self.define_solution_space(p)

         address = self._read_address(p)
         metfile, modlen = p.gen.metfile, 3*len(address)
         
         ### get gmat, Cd and d values and send segments to other processors  
         PETSc.Sys.Print('retrieving gmat, cdi, and dvec vals...')
         gmatvals, gmatij, cdiv_temp, dvec_temp = get_gcd_vals(metfile,len(p.scenes),
                  p.gen.doline,p.gen.cols,address,p.gen.gmat_rows)

         ### estimate posterior model pixel x pixel (for PETSc initial guess)
         PETSc.Sys.Print('estimating posterior model...')
         mtilde_temp = pixelxpixelsolver(metfile,len(p.scenes),
                  p.gen.doline,p.gen.cols,address)
         gmatij = gmatij.reshape(-1,3*p.gen.gmat_rows)
         if gmatij.dtype != np.int32:  gmatij = gmatij.astype(np.int32)
         if p.gen.cmlambda < 1.e-12:  del address
        
         ### initialize arrays to store start and end row indices for all procs 
         rowblocks_gmat = np.empty([2,size],dtype=np.int32)
         rowblocks_cdiv = np.empty([2,size],dtype=np.int32)
         rowblocks_dvec = np.empty([2,size],dtype=np.int32)
         rowblocks_mtil = np.empty([2,size],dtype=np.int32)
      
         tag = 1234
         for i in np.arange(1,size):
            comm.Send([np.array([modlen,p.gen.gmat_rows]),MPI.INT], dest=i, tag=tag); tag += 1
      else:
         tag, pararray = 1234 + (rank-1), np.array([-1,-1])
         comm.Recv([pararray,MPI.INT], source=0, tag=tag)
         modlen, p.gen.gmat_rows = pararray[0], pararray[1]
         if modlen < 0 or p.gen.gmat_rows < 0:
            print('\npararray receive command failed in rank '+str(rank))
            PETSc._finalize; MPI.Finalize; sys.exit()

      PETSc.Sys.Print('building gmat, cdiv, and dvec arrays...')
      ### initialize...nnz --> [DIAGONAL,OFF-DIAGONAL] (see PETSc documentation)
      gmat = PETSc.Mat().createAIJ([p.gen.gmat_rows,modlen],nnz=[3,3],comm=comm)
      cdiv = PETSc.Mat().createAIJ([p.gen.gmat_rows,p.gen.gmat_rows],nnz=[1,1],comm=comm)
      dvec = PETSc.Vec().createMPI(p.gen.gmat_rows,comm=comm)
      mtilde_pre = PETSc.Vec().createMPI(modlen,comm=comm)
      ### get the block of rows owned by this processor 
      sr_gmat, er_gmat = gmat.getOwnershipRange()
      sr_cdiv, er_cdiv = cdiv.getOwnershipRange()
      sr_dvec, er_dvec = dvec.getOwnershipRange()
      sr_mtil, er_mtil = mtilde_pre.getOwnershipRange()

      ### send start/end rows to root proc and receive gmat, cdiv, and dvec entries 
      base_tag, base_stag = 777, 999
      if rank == 0:
         rowblocks_gmat[0,0], rowblocks_gmat[1,0] = sr_gmat, er_gmat
         rowblocks_cdiv[0,0], rowblocks_cdiv[1,0] = sr_cdiv, er_cdiv
         rowblocks_dvec[0,0], rowblocks_dvec[1,0] = sr_dvec, er_dvec
         rowblocks_mtil[0.0], rowblocks_mtil[1,0] = sr_mtil, er_mtil
         tag,stag = base_tag, base_stag
         recsimple = np.array([-1,-1])
         for i in np.arange(1,size):
            comm.Recv([recsimple,MPI.INT],source=i,tag=tag); tag += 1
            rowblocks_gmat[0,i], rowblocks_gmat[1,i] = recsimple[0], recsimple[1]
            comm.Recv([recsimple,MPI.INT],source=i,tag=tag); tag += 1
            rowblocks_cdiv[0,i], rowblocks_cdiv[1,i] = recsimple[0], recsimple[1]
            comm.Recv([recsimple,MPI.INT],source=i,tag=tag); tag += 1
            rowblocks_dvec[0,i], rowblocks_dvec[1,i] = recsimple[0], recsimple[1]
            comm.Recv([recsimple,MPI.INT],source=i,tag=tag); tag += 1
            rowblocks_mtil[0,i], rowblocks_mtil[1,i] = recsimple[0], recsimple[1]

            svec = gmatvals[3*rowblocks_gmat[0,i]:3*rowblocks_gmat[1,i]]  
            comm.Send([svec, MPI.FLOAT], dest=i, tag=stag); stag += 1

            svec = gmatij[:,3*rowblocks_gmat[0,i]:3*rowblocks_gmat[1,i]].flatten()
            comm.Send([svec, MPI.INT], dest=i, tag=stag); stag += 1

            svec = cdiv_temp[rowblocks_cdiv[0,i]:rowblocks_cdiv[1,i]]
            comm.Send([svec, MPI.FLOAT], dest=i, tag=stag); stag += 1

            svec = dvec_temp[rowblocks_dvec[0,i]:rowblocks_dvec[1,i]]
            comm.Send([svec, MPI.FLOAT], dest=i, tag=stag); stag += 1

            svec = mtilde_temp[rowblocks_mtil[0,i]:rowblocks_mtil[1,i]]
            comm.Send([svec, MPI.FLOAT], dest=i, tag=stag); stag += 1

         if size > 1: del svec

         gmatvals = gmatvals[3*rowblocks_gmat[0,0]:3*(rowblocks_gmat[1,0])]
         gmatij = gmatij[:,3*rowblocks_gmat[0,0]:3*(rowblocks_gmat[1,0])]
         cdiv_temp = cdiv_temp[rowblocks_cdiv[0,0]:rowblocks_cdiv[1,0]]
         dvec_temp = dvec_temp[rowblocks_dvec[0,0]:rowblocks_dvec[1,0]]
         mtilde_temp = mtilde_temp[rowblocks_mtil[0,0]:rowblocks_mtil[1,0]]

      else:
         tag = base_tag + (rank-1)*4
         comm.Send([np.array([sr_gmat, er_gmat]),MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([np.array([sr_cdiv, er_cdiv]),MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([np.array([sr_dvec, er_dvec]),MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([np.array([sr_mtil, er_mtil]),MPI.INT], dest=0, tag=tag); tag += 1

         stag = base_stag + (rank-1)*5

         gmatvals = np.empty(3*(er_gmat-sr_gmat), dtype=np.float32)
         gmatij = np.empty(2*3*(er_gmat-sr_gmat), dtype=np.int32)
         cdiv_temp = np.empty((er_cdiv-sr_cdiv), dtype=np.float32)
         dvec_temp = np.empty((er_dvec-sr_dvec), dtype=np.float32)
         mtilde_temp = np.empty((er_mtil-sr_mtil), dtype=np.float32)

         comm.Recv([gmatvals,MPI.FLOAT],source=0, tag=stag); stag += 1
         comm.Recv([gmatij,MPI.INT],source=0, tag=stag); stag += 1
         comm.Recv([cdiv_temp,MPI.FLOAT],source=0, tag=stag); stag += 1
         comm.Recv([dvec_temp,MPI.FLOAT],source=0, tag=stag); stag += 1
         comm.Recv([mtilde_temp,MPI.FLOAT],source=0, tag=stag); stag += 1
         gmatij = gmatij.reshape(2,-1)

      PETSc.Sys.Print('loading gmat, cdiv, and dvec sparse matrices...')
      for i in np.arange(len(dvec_temp)):
         threei = 3*i
         gmat.setValue(gmatij[0,threei], gmatij[1,threei], gmatvals[threei])
         gmat.setValue(gmatij[0,threei], gmatij[1,threei+1], gmatvals[threei+1])
         gmat.setValue(gmatij[0,threei], gmatij[1,threei+2], gmatvals[threei+2])
         cdiv.setValue(gmatij[0,threei], gmatij[0,threei], cdiv_temp[i])
         dvec.setValue(gmatij[0,threei], dvec_temp[i])
      comm.Barrier()
      dvec.assemblyBegin(); dvec.assemblyEnd()  # MAT_FINAL_ASSEMBLY is implied 
      gmat.assemblyBegin(); gmat.assemblyEnd()
      cdiv.assemblyBegin(); cdiv.assemblyEnd()
      comm.Barrier()
      for i in np.arange(len(mtilde_temp)):
         mtilde_pre.setValue(i+sr_mtil, mtilde_temp[i])
      comm.Barrier()
      mtilde_pre.assemblyBegin(); mtilde_pre.assemblyEnd()
      comm.Barrier()
      del gmatvals, gmatij, cdiv_temp, dvec_temp, mtilde_temp

      ### build A and gtilde (=dtilde)

      omega = gmat.transposeMatMult(cdiv)
      cdiv.destroy()

      gtg = omega.matMult(gmat)
      mtilde,gtilde = gtg.getVecs() 
      omega.mult(dvec,gtilde)
      mtilde_pre.copy(mtilde)
      dvec.destroy(); omega.destroy(); gmat.destroy(); mtilde_pre.destroy()

      if p.gen.cmlambda > 1.e-12:

         PETSc.Sys.Print('building Laplacian...')

         if rank == 0:
            try:
               aa = len(address)
               if 3*aa != modlen: 
                  address = self._read_address(p)
            except:
               address = self._read_address(p)
            fmatval, fmatij = laplace_builder(p.gen.cols,address)

            fmatval = fmatval[:np.argmin(fmatval)]
            fmatij = fmatij.reshape(2,-1)[:,:len(fmatval)]
            if fmatij.dtype != np.int32: fmatij = fmatij.astype(np.int32)
            del address

            rowblocks_fmat = np.empty([2,size],dtype=np.int32)

         fmat = PETSc.Mat().createAIJ([modlen, modlen],nnz=[5,5],comm=comm) 

         ### get the row addresses for this processor 
         sr_gtg, er_gtg = gtg.getOwnershipRange()
         sr_fmat, er_fmat = fmat.getOwnershipRange()

         base_tag, base_stag = 77, 444
         if rank == 0:
            rowblocks_fmat[0,0], rowblocks_fmat[1,0] = sr_fmat, er_fmat
            tag, stag = base_tag, base_stag
            recarray = np.array([-1,-1])
            for i in np.arange(1,size):
               comm.Recv([recarray,MPI.INT],source=i,tag=tag); tag += 1
               rowblocks_fmat[:,i] = recarray
               # get indexes for fmatval and fmatij
               left  = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[0,i]))
               right = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[1,i]))

               svec = fmatval[left:right]
               lensvec = np.array([len(svec)],dtype=np.int32)

               comm.Send([lensvec,MPI.INT], dest=i, tag=stag); stag += 1
               comm.Send([svec,MPI.FLOAT], dest=i, tag=stag); stag += 1                

               svec = fmatij[:,left:right].flatten()
               comm.Send([svec,MPI.INT], dest=i, tag=stag); stag += 1

            if size > 1: del svec
            left  = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[0,0]))
            right = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[1,0]))
            fmatval = fmatval[left:right]
            fmatij = fmatij[:,left:right]

         else:
            tag, stag = base_tag + (rank-1), base_stag + (rank-1)*3
            lenrec = np.array([-1],dtype=np.int32)
            comm.Send([np.array([sr_fmat, er_fmat]),MPI.INT], dest=0, tag=tag); tag += 1
            comm.Recv([lenrec,MPI.INT], source=0, tag=stag); stag += 1
            lenrec = lenrec[0]

            fmatval = np.empty(lenrec,dtype=np.float32)
            fmatij  = np.empty(2*lenrec,dtype=np.int32)

            comm.Recv([fmatval,MPI.FLOAT], source=0, tag=stag); stag += 1
            comm.Recv([fmatij,MPI.INT], source=0, tag=stag); stag += 1 
            fmatij = fmatij.reshape(2,-1)


         PETSc.Sys.Print('loading sparse Laplacian...')
         for i in np.arange(len(fmatval)):
            fmat.setValue(fmatij[0,i],fmatij[1,i],fmatval[i])
         comm.Barrier()
         fmat.assemblyBegin(); fmat.assemblyEnd()
         comm.Barrier()
         del fmatval, fmatij

         PETSc.Sys.Print('loading cli...')
         cliv = PETSc.Vec().createMPI(modlen,comm=comm)
         gtg.getDiagonal(cliv)
         cli = gtg.duplicate(copy=False)
         cli.setDiagonal(cliv)
         cliv.destroy()

         PETSc.Sys.Print('building C_m...')
         rhs = fmat.transposeMatMult(cli)
         rhs = rhs.matMult(fmat)
         cli.destroy(); fmat.destroy()

         PETSc.Sys.Print('building A...')
         gtg.axpy(p.gen.cmlambda,rhs)
         rhs.destroy()

      #if p.gen.outputs['Output diag(G^T*W*G)^-1)']:
         #PETSc.Sys.Print('Output diag(G^T*W*G)^-1) option not implemented...run vector_disp for estimate')
      '''
         ### this part of the code works, but the inverse operator is not implemented due to expense 
         PETSc.Sys.Print('writing diagonal of posterior model covariance matrix')
         cliv = PETSc.Vec().createMPI(modlen,comm=comm)
         gtg.getDiagonal(cliv)
         
         base_tag, base_stag = 2222, 3333
         sr_mt, er_mt = cliv.getOwnershipRange()   
         if rank == 0:
            clivnump = np.empty(modlen, dtype=np.float32)
            clivnump[sr_mt:er_mt] = cliv[...]; cliv.destroy()

            rowarr = np.array([-1,-1], dtype=np.int32)  
            tag = base_tag
            for i in np.arange(1,size):
               comm.Recv([rowarr,MPI.INT], source=i, tag=tag); tag += 1
               comm.Recv([clivnump[rowarr[0]:rowarr[1]], MPI.FLOAT], source=i, tag=tag); tag += 1

            self._write_post_model(model=clivnump,p=p,form=1)
         else:
            rowarr, tag = np.array([sr_mt, er_mt], dtype=np.int32), base_tag + (rank-1)*2
            clivnump = cliv[...]; cliv.destroy()

            if clivnump.dtype != np.float32:  clivnump = clivnump.astype(np.float32)
            comm.Send([rowarr,MPI.INT], dest=0, tag=tag); tag += 1
            comm.Send([clivnump, MPI.FLOAT], dest=0, tag=tag); tag += 1
      '''

      PETSc.Sys.Print('inverting using PETSc lsqr...') 
      ksp = PETSc.KSP().create(comm)
      ksp.setType('lsqr')
      ksp.pc.setType('none')
      ksp.setOperators(gtg) 
      self._get_ksp_tols(ksp=ksp,p=p)
      ksp.setTolerances(rtol=self.rtol,atol=self.atol,divtol=self.dtol,max_it=self.max_it)
      t0 = time.time()
      ksp.solve(gtilde,mtilde)
      tf = time.time()

      gtg.destroy(); gtilde.destroy()

      PETSc.Sys.Print('  lsqr took %10.2f seconds' % (tf-t0))
      PETSc.Sys.Print('  Converged in %d iterations ' % (ksp.getIterationNumber()))
      PETSc.Sys.Print('  Tolerances: %e %e %d %d' % (ksp.getTolerances()))
      PETSc.Sys.Print('  Convergance code (< 0 --> diverged): %d\n' % (ksp.getConvergedReason()))

      ksp.destroy()

      base_tag, base_stag = 2222, 3333
      sr_mt, er_mt = mtilde.getOwnershipRange() 
      if rank == 0:    
         mtildenump = np.empty(modlen, dtype=np.float32)
         mtildenump[sr_mt:er_mt] = mtilde[...]; mtilde.destroy()
         mtildenump *= p.gen.timefac * p.gen.distfac * p.gen.r2dsp

         rowarr = np.array([-1,-1], dtype=np.int32)
         tag = base_tag
         for i in np.arange(1,size):
            comm.Recv([rowarr,MPI.INT], source=i, tag=tag); tag += 1
            comm.Recv([mtildenump[rowarr[0]:rowarr[1]], MPI.FLOAT], source=i, tag=tag); tag += 1
         self._write_post_model(model=mtildenump,p=p)

         self._cleanup(p=p)
      else:
         rowarr, tag = np.array([sr_mt, er_mt], dtype=np.int32), base_tag + (rank-1)*2
         mtildenump = mtilde[...]; mtilde.destroy()
         mtildenump *= p.gen.timefac * p.gen.distfac * p.gen.r2dsp

         if mtildenump.dtype != np.float32:  mtildenump = mtildenump.astype(np.float32) 
         comm.Send([rowarr,MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([mtildenump, MPI.FLOAT], dest=0, tag=tag); tag += 1

      PETSc._finalize
      MPI.Finalize
###-----------------------------------------------------------------------------

 
###=============================================================================

### private functions 

###=============================================================================
   def _get_ksp_tols(self,ksp,p):
      """
      Tolerances for PETSc KSP solver
         rtol   : relative convergence tolerance [1e-5]
         atol   : absolute convergence tolerance [1e-50]
         dtol   : divergence tolerance
         max_it : max number of interations 
      #if rank==0: print('ksp.getTolerances(): ',ksp.getTolerances())
      """
      defaults = ksp.getTolerances()
      self.rtol, self.atol = p.gen.petsc_rtol, p.gen.petsc_atol
      self.dtol, self.max_it = p.gen.petsc_dtol, p.gen.petsc_max_it
      if not self.rtol:   self.rtol = defaults[0]
      if not self.atol:   self.atol = defaults[1]
      if not self.dtol:   self.dtol = defaults[2]
      if not self.max_it: self.max_it = defaults[3]

###-----------------------------------------------------------------------------
   def _write_metfile(self,p,filename='.r3d.metfile',delim=':'):
      fid = open(filename,'w')
      lp = len(p.scenes)
      for i in np.arange(lp):
         wstr  = p.scenes[i].unwfile + delim + p.scenes[i].corfile + delim
         wstr += p.scenes[i].losfile + delim + str(p.scenes[i].cols) + delim
         wstr += str(p.scenes[i].phznull) + delim + str(p.scenes[i].cornull) + delim
         wstr += str(p.scenes[i].losnull) + delim + str(p.scenes[i].offlat) + delim
         wstr += str(p.scenes[i].sline) + delim + str(p.scenes[i].eline) + delim
         wstr += str(p.scenes[i].bcol) + delim + str(p.scenes[i].scol) + delim
         wstr += str(p.scenes[i].dlooks) + delim + str(p.scenes[i].dsp2r) + '\n'
         fid.write(wstr)
      fid.close()

###-----------------------------------------------------------------------------
   def _write_post_model(self,model,p,form=0):
      """
      Writes the posterior model to separate east, north, up files

      form = 0   for mtilde (velocity)
      form = 1   for posterior model covariance diagonal components (*.cov.*) 
      """
      outname = p.gen.outfldr + p.gen.outpref
      if form == 1:
         outname += '.cov' 
      address  = self._read_address(p)
      outeast  = np.empty([p.gen.doline,p.gen.cols],dtype=np.float32)
      outnorth = np.empty([p.gen.doline,p.gen.cols],dtype=np.float32)
      outup    = np.empty([p.gen.doline,p.gen.cols],dtype=np.float32)
      addind, lenadd = 0, len(address)
      for i in np.arange(p.gen.doline):
         for j in np.arange(p.gen.cols):
            if addind < lenadd and address[addind] == i*p.gen.cols+j:
               outeast[i,j] = model[3*addind]
               outnorth[i,j] = model[3*addind+1]
               outup[i,j] = model[3*addind+2]
               addind += 1
            else:
               outeast[i,j] = p.gen.outnull
               outnorth[i,j] = p.gen.outnull
               outup[i,j] = p.gen.outnull
      fid = open(outname+'.east','w')
      outeast.flatten().tofile(fid)
      fid.close()
      fid = open(outname+'.north','w')
      outnorth.flatten().tofile(fid)
      fid.close()
      fid = open(outname+'.up','w')
      outup.flatten().tofile(fid)
      fid.close() 

###-----------------------------------------------------------------------------
   def _write_numscenes(self,numscenes,p):
      fid = open(p.gen.numfile,'w')
      numscenes.flatten().tofile(fid)
      fid.close()
   def _read_numscenes(self,p):
      fid = open(p.gen.numfile,'r')
      numscenes = np.fromfile(fid,dtype=np.float32)
      fid.close()
      return numscenes

###-----------------------------------------------------------------------------
   def _write_address(self,address,p):
      fid = open(p.gen.addfile,'w')
      address.astype(np.int32).tofile(fid)
      fid.close()
   def _read_address(self,p):
      fid = open(p.gen.addfile,'r')
      address = np.fromfile(fid,dtype=np.int32)
      fid.close()
      return address

###-----------------------------------------------------------------------------
   def _cleanup(self,p):
      os.remove(p.gen.addfile)  # remove address file
      os.remove(p.gen.metfile)  # remove metfile
      os.remove('.gmat_rows')

###-----------------------------------------------------------------------------
   def _open_datafiles(self,p):
      """
      Open files for reading
      """
      fidunw, fidcor, fidlos = [], [], []
      for i in np.arange(len(p.scenes)):
         fidunw.append(open(p.scenes[i].unwfile,'r'))
         fidcor.append(open(p.scenes[i].corfile,'r'))
         fidlos.append(open(p.scenes[i].losfile,'r'))

         if p.scenes[i].offlat < 0:  # set cursor 
            fidunw[i].seek(-p.scenes[i].offlat*p.scenes[i].cols*4,0)
            fidcor[i].seek(-p.scenes[i].offlat*p.scenes[i].cols*4,0)
            fidlos[i].seek(-p.scenes[i].offlat*p.scenes[i].cols*12,0)
      return fidunw,fidcor,fidlos

###-----------------------------------------------------------------------------
   def _read_data_lines(self,unw,cor,los,dum,fidunw,fidcor,fidlos,p,i):
      cnt = 0 
      for m in np.arange(len(p.scenes)):
         if i >= p.scenes[m].offlat and i <= p.scenes[m].eline:
            cnt += 1
            cols, scol, bcol = p.scenes[m].cols, p.scenes[m].scol, p.scenes[m].bcol
 
            dum[:] = p.scenes[m].phznull 
            dum[bcol:bcol+cols] = np.fromfile(fidunw[m],dtype=np.float32,count=cols)  
            unw[m,:] = dum[scol:p.gen.cols+scol]

            if p.scenes[m].phznull != p.scenes[m].cornull:  dum[:] = p.scenes[m].cornull
            dum[bcol:bcol+cols] = np.fromfile(fidcor[m],dtype=np.float32,count=cols)  
            cor[m,:] = dum[scol:p.gen.cols+scol]

            if p.scenes[m].cornull != p.scenes[m].losnull: dum[:] = p.scenes[m].losnull
            dum[3*bcol:3*(bcol+cols)] = np.fromfile(fidlos[m],dtype=np.float32,count=3*cols)
            los[m,:] = dum[3*scol:3*(p.gen.cols+scol)]
         else:
            unw[m,:] = p.scenes[m].phznull
      return unw,cor,los,cnt

###-----------------------------------------------------------------------------
   def _ffwd_data(self,fidunw,fidcor,fidlos,p,i):
      for m in np.arange(len(p.scenes)):
         if i >= p.scenes[m].offlat and i < p.scenes[m].eline:
            fidunw[m].seek(p.scenes[m].cols*4,1)
            fidcor[m].seek(p.scenes[m].cols*4,1)
            fidlos[m].seek(p.scenes[m].cols*12,1)  

###-----------------------------------------------------------------------------
   def _write_output(self,data,p,ext):
      outname = p.gen.outfldr + p.gen.outpref + ext
      fid = open(outname,'w')
      data.flatten().tofile(fid)
      fid.close()

###-----------------------------------------------------------------------------
   def _interograte_data(self,unw,cor,los,p,m,nulltol=1.e-8):
      """
      Returns True if no null values are found, no NaNs are present, and all values are finite
      """
      dumb  = np.any(np.abs(los-p.scenes[m].losnull) < nulltol)
      dumb += np.abs(unw - p.scenes[m].phznull) < nulltol
      dumb += np.abs(cor - p.scenes[m].cornull) < nulltol
      dumb += np.any(np.isnan(los))
      dumb += np.isnan(unw)
      dumb += np.isnan(cor)
      dumb += np.any(np.isinf(los))
      dumb += np.isinf(unw)
      dumb += np.isinf(cor)
      dumb += np.abs(np.linalg.norm(los) - 1.) > 1.e-3
      return -dumb  

###-----------------------------------------------------------------------------
   def _initialize_input_vecs(self,p):
      lp = len(p.scenes)
      unw = np.empty([lp,p.gen.cols],dtype=np.float32)
      cor = np.empty([lp,p.gen.cols],dtype=np.float32)
      los = np.empty([lp,3*p.gen.cols],dtype=np.float32)
      dum = np.empty(10*p.gen.cols,dtype=np.float32)
      return unw,cor,los,dum


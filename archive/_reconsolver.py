# -*- coding: utf-8 -*-
"""

Module to read data and load the design (G) and covariance (C_phi) matrices 

"""

import datetime
import sys,os,time
import numpy as np
from scipy.sparse import dok_matrix
import _recontools, _reconparser

__all__ = ['Solver']

###===========================================================

###===========================================================


class Solver():
   """
   Solver Class

   Contains all of the solver functions.  This is the main engine for the fluye code.

   Usage:  Solver(parser_class)  
   """

   def __init__(self,p):
      sys.stdout.write('Instantiating Solver class\n')
      self.note = 'Created on', datetime.datetime.now()

   def define_solution_space(self,p):
      from _reconutils import defsolspace
      from scipy.linalg import svdvals
      sys.stdout.write('Defining the solution space\n')
      metfile, delim = '.fluye.metfile', ':'
      self._write_metfile(p,filename=metfile,delim=delim)
      lp = len(p.scenes)

      unw,cor,los,dum = self._initialize_input_vecs(p=p)
      fidunw,fidcor,fidlos = self._open_datafiles(p=p)

      ###  loop through the image to find the addresses of the non-Null space pixels
      address  = np.empty(p.gen.doline*p.gen.cols,dtype=np.int32)
      nonulout = np.empty([p.gen.doline,p.gen.cols],dtype=np.float32)
      rankout  = np.empty([p.gen.doline,p.gen.cols],dtype=np.float32)
      losarray = np.empty([lp,3])

      addind = 0   # index for address vector
      p.gen.gmat_rows = 0
      colvec, scvec = np.arange(p.gen.cols), np.arange(lp)

      for i in np.arange(p.gen.doline):
         percent = 100*i//p.gen.doline
         if percent % 2 == 0: 
            sys.stdout.write("%4d%%\r" % percent)
            sys.stdout.flush()
            
         unw,cor,los,cnt = self._read_data_lines(unw,cor,los,dum,fidunw,fidcor,fidlos,p,i) 
         if cnt >= 3:
            for j in colvec:
               nonul, jind = 0, 3*j
               for m in scvec:
                  if self._interograte_data(unw=unw[m,j],cor=cor[m,j],los=los[m,jind:jind+3],p=p,m=m):
                     losarray[nonul,:] = los[m,jind:jind+3]
                     nonul += 1
               if nonul > 2:
                  nonulout[i,j] = nonul
                  # determine linear independence of los vectors 
                  svdv = svdvals(losarray[:nonul],check_finite=False)
                  losrank = np.sum( svdv > 1.e-7 )
                  rankout[i,j] = losrank.astype(np.float32)
                  if losrank > 2:
                     address[addind] = i*p.gen.cols+j
                     addind += 1
                     p.gen.gmat_rows += nonul
               else:
                  nonulout[i,j], rankout[i,j] = p.gen.outnull, p.gen.outnull
         else:
            nonulout[i], rankout[i] = p.gen.outnull, p.gen.outnull
      del colvec
      percent = 100
      sys.stdout.write('%4d%%\n\n' % percent)
      sys.stdout.write('Writing addresses of nonnull values\n')
      
      address = address[:addind]
      self._write_address(address=address,p=p)

      fid = open('.gmat_rows','w')
      fid.write(str(p.gen.gmat_rows))
      fid.close()

      ### close data files 
      for m in np.arange(lp):
         fidunw[m].close(); fidcor[m].close(); fidlos[m].close()
      ### write out num and rank files if requested 
      if p.gen.outputs['Output number of scenes']:
         self._write_output(data=nonulout,p=p,ext='.num')
      if p.gen.outputs['Output rank(G) at each pixel']:
         self._write_output(data=rankout,p=p,ext='.rnk')


      
   def BS(self,p):
      """

      Basic solver that inverts the entire scene using a single patch.
      The inversion is carried out with scipy.sparse.linalg.lsqr

      """

      sys.stdout.write("\nRunning BS:  Basic solver using a single patch\n")
      from scipy.sparse import dok_matrix, dia_matrix, bsr_matrix, coo_matrix, identity
      from scipy.sparse.linalg import lsqr

      cols = p.gen.cols
      address = self._read_address(p)
      unw,cor,los,dum = self._initialize_input_vecs(p=p)
      fidunw,fidcor,fidlos = self._open_datafiles(p) 

      try:
         p.gen.gmat_rows
      except:
         fid = open('.gmat_rows','r')
         p.gen.gmat_rows = np.int32(''.join(fid.readlines()))
         fid.close()


      sys.stdout.write('building sparse matrices\n')
      addind, gmatrow, ott = 0, 0, np.array([0,1,2])
      colvec, scvec = np.arange(cols), np.arange(len(p.scenes))
      gmatvals = np.empty(p.gen.gmat_rows*3,dtype=np.float64)
      gmatij = np.empty([2,p.gen.gmat_rows*3],dtype=np.int32)
      cdiv = np.empty(p.gen.gmat_rows,dtype=np.float64)
      dvec = np.empty(p.gen.gmat_rows,dtype=np.float64)

      for i in np.arange(p.gen.doline):
         percent = 100*i//p.gen.doline
         if percent % 2 == 0: 
            sys.stdout.write("%4d%%\r" % percent)
            sys.stdout.flush()

         if address[addind] >= i*cols or address[addind] < (i+1)*cols: # check if there are good vals in row
            unw,cor,los,cnt = self._read_data_lines(unw,cor,los,dum,fidunw,fidcor,fidlos,p,i)
            for j in colvec:
               if address[addind] == i*cols+j:
                  jind = 3*j
                  for m in scvec:
                     if self._interograte_data(unw=unw[m,j],cor=cor[m,j],los=los[m,jind:jind+3],p=p,m=m):
                        gind = 3*gmatrow
                        corsq = cor[m,j]**2
                        gmatij[0,gind:gind+3] = gmatrow
                        gmatij[1,gind:gind+3] = 3*addind + ott
                        gmatvals[gind:gind+3] = los[m,jind:jind+3]
                        dvec[gmatrow] = unw[m,j] * p.scenes[m].dsp2r
                        cdiv[gmatrow] = p.scenes[m].dlooks*corsq/(1.-corsq)
                        gmatrow += 1
                  addind += 1
                  if addind >= len(address):  addind -= 1
         else:
            self._ffwd_data(fidunw=fidunw,fidcor=fidcor,fidlos=fidlos,p=p,i=i)
      percent = 100
      sys.stdout.write("%4d%%\n\n" % percent)
      
      # stow data files and release memory       
      del unw,cor,los,dum,colvec,scvec
      for m in np.arange(len(p.scenes)):
         fidunw[m].close(); fidcor[m].close(); fidlos[m].close()

      inds = np.arange(p.gen.gmat_rows)

      gmat = bsr_matrix((gmatvals,gmatij),dtype=np.float64)
      cdi  = coo_matrix((cdiv,(inds,inds)),dtype=np.float64).todia()
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
      fmat = self._build_laplace_scipy_sparse(a=address,p=p)
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
            mtilde = lsqr(A,dtilde,damp=p.gen.damp,conlim=1.e8,iter_lim=2e5,show=False)  
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

###==============================================================================

   def _write_metfile(self,p,filename='.fluye.metfile',delim=':'):
      fid = open(filename,'w')
      lp = len(p.scenes)
      for i in np.arange(lp):
         wstr  = p.scenes[i].unwfile + delim + p.scenes[i].corfile + delim
         wstr += p.scenes[i].losfile + delim + str(p.scenes[i].cols) + delim
         wstr += str(p.scenes[i].phznull) + delim + str(p.scenes[i].cornull) + delim
         wstr += str(p.scenes[i].losnull) + delim + str(p.scenes[i].offlat) + delim
         wstr += str(p.scenes[i].sline) + delim + str(p.scenes[i].eline) + delim
         wstr += str(p.scenes[i].bcol) + delim + str(p.scenes[i].scol) + '\n'
         fid.write(wstr)
      fid.close()

   def _write_post_model(self,model,p):
      """
      Writes the posterior model to separate east, north, up files
      """
      outname = p.gen.outfldr + p.gen.outpref 
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
       

   def _build_laplace_scipy_sparse(self,a,p):
      """
      Builds a Laplacian for the model vector 

      Parameters
      ----------
      a     :     ndarray
                  Array of addresses for the model components 
      p     :     Parser class
                  Contains info parsed from the command file

      """
      laplace = dok_matrix((3*len(a),3*len(a)),dtype=np.float64)
      cols = p.gen.cols
      for i in np.arange(len(a)):
         row, diag = 3*i, 0.
         left, right, above, below = a[i]-1, a[i]+1, a[i]-cols, a[i]+cols
         # left must be positive and i > 0 for there to be a valid pixel to the left
         if left >= 0 and i > 0:
            # the way address are loaded means that if the pixel to the left is valid, it 
            # will be at i-1. (Address is row+col.) For the pixel to be to the left of the 
            # current pixel, the current pixel can't be on the left image boundary
            if a[i-1] == left and a[i] % cols != 0:
               j = 3*(i-1)
               laplace[row,j], laplace[row+1,j+1], laplace[row+2,j+2] = -1., -1., -1.
               diag += 1.
         # basically the same arguments as for the left pixel
         if right < cols and i < len(a)-1:
            if a[i+1] == right and a[i+1] % cols != 0:
               j = 3*(i+1)
               laplace[row,j], laplace[row+1,j+1], laplace[row+2,j+2] = -1., -1., -1.
               diag += 1.
         # Address must put the current pixel in at least the 2nd row 
         if above >= 0:
            sl = np.max([i-(cols+1),0])
            # find out if above pixel is in address list...no need to search further 
            # back than i-cols because of how address vector was loaded
            if np.any(a[sl:i]==above):
               j = 3*( sl + np.argmin(np.abs(a[sl:i]-above)))
               laplace[row,j], laplace[row+1,j+1], laplace[row+2,j+2] = -1., -1., -1.
               diag += 1.
         if below < p.gen.doline*cols:  # must not be greater than the highest address
            sl = np.min([i+cols+1,len(a)])
            if np.any(a[i:sl]==below):
               j = 3*( i + np.argmin(np.abs(a[i:sl]-below)))
               laplace[row,j], laplace[row+1,j+1], laplace[row+2,j+2] = -1., -1., -1.
               diag += 1.
         laplace[row,row], laplace[row+1,row+1], laplace[row+2,row+2] = diag, diag, diag
      return laplace
 
      

   def _write_numscenes(self,numscenes,p):
      fid = open(p.gen.numfile,'w')
      numscenes.flatten().tofile(fid)
      fid.close()
   def _read_numscenes(self,p):
      fid = open(p.gen.numfile,'r')
      numscenes = np.fromfile(fid,dtype=np.float32)
      fid.close()
      return numscenes

   def _write_address(self,address,p):
      fid = open(p.gen.addfile,'w')
      address.astype(np.int32).tofile(fid)
      fid.close()
   def _read_address(self,p):
      fid = open(p.gen.addfile,'r')
      address = np.fromfile(fid,dtype=np.int32)
      fid.close()
      return address

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


   def _ffwd_data(self,fidunw,fidcor,fidlos,p,i):
      for m in np.arange(len(p.scenes)):
         if i >= p.scenes[m].offlat and i < p.scenes[m].eline:
            fidunw[m].seek(p.scenes[m].cols*4,1)
            fidcor[m].seek(p.scenes[m].cols*4,1)
            fidlos[m].seek(p.scenes[m].cols*12,1)  



   def _write_output(self,data,p,ext):
      outname = p.gen.outfldr + p.gen.outpref + ext
      fid = open(outname,'w')
      data.flatten().tofile(fid)
      fid.close()

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

   def _initialize_input_vecs(self,p):
      lp = len(p.scenes)
      unw = np.empty([lp,p.gen.cols],dtype=np.float32)
      cor = np.empty([lp,p.gen.cols],dtype=np.float32)
      los = np.empty([lp,3*p.gen.cols],dtype=np.float32)
      dum = np.empty(10*p.gen.cols,dtype=np.float32)
      return unw,cor,los,dum


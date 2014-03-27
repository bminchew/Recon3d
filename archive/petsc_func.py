
   def PPS(self,p):
      """
      Parallel solver implementing a PETSc KSP linear routine over the entire scene as a 
      single patch (PPS = Parallel PETSc Single)
      """
      from _reconutils import get_gcd_vals, laplace_builder
      try:
         import mpi4py
         from mpi4py import MPI
      except:
         print 'PPS requires that mpi4py be installed and in the Python path'
         sys.exit()
      try:
         import petsc4py
         petsc4py.init([],comm=MPI.COMM_WORLD)
         from petsc4py import PETSc
      except:
         print 'PPS requires that petsc4py be installed and in the Python path'
         sys.exit()
 
      PETSc.Sys.Print('Running PPS:  Parallel PETSc (linear) Solver on a single patch')

      comm=PETSc.COMM_WORLD
      size, rank = comm.Get_size(), comm.Get_rank()
      cols = p.gen.cols
      try:
         p.gen.gmat_rows
      except:
         fid = open('.gmat_rows','r')
         p.gen.gmat_rows = np.int32(''.join(fid.readlines()))
         fid.close()

      if rank == 0:
         address = self._read_address(p)
         metfile, modlen = p.gen.metfile, 3*len(address)
         PETSc.Sys.Print('retrieving gmat, cdi, and dvec vals')
         
         ### get gmat, Cd and d values and send segments to other processors  
         gmatvals, gmatij, cdiv_temp, dvec_temp = get_gcd_vals(metfile,len(p.scenes),
                  p.gen.doline,p.gen.cols,address,p.gen.gmat_rows)
         gmatij = gmatij.reshape(-1,3*p.gen.gmat_rows)
         if gmatij.dtype != np.int32:  gmatij = gmatij.astype(np.int32)
         if p.gen.cmlambda < 1.e-12:  del address
        
         ### initialize arrays to store start and end row indices for all procs 
         rowblocks_gmat = np.empty([2,size],dtype=np.int32)
         rowblocks_cdiv = np.empty([2,size],dtype=np.int32)
         rowblocks_dvec = np.empty([2,size],dtype=np.int32)

      modlen = comm.bcast(modlen,root=0)

      PETSc.Sys.Print('building sparse matrices')
      ### initialize 
      gmat = PETSc.Mat().createAIJ([p.gen.gmat_rows,modlen],nnz=p.gen.gmat_rows*3,comm=comm)
      cdiv = PETSc.Mat().createAIJ([p.gen.gmat_rows,p.gen.gmat_rows],nnz=p.gen.gmat_rows,comm=comm)
      dvec = PETSc.Vec().createMPI(p.gen.gmat_rows,comm=comm)
         
      ### get the block of rows owned by this processor 
      sr_gmat, er_gmat = gmat.getOwnershipRange()
      sr_cdiv, er_cdiv = cdiv.getOwnershipRange()
      sr_dvec, er_dvec = dvec.getOwnershipRange()

      ### send start/end rows to root proc and receive gmat, cdiv, and dvec entries 
      base_tag, base_stag = 777, 999
      if rank == 0:
         rowblocks_gmat[0,0], rowblocks_gmat[1,0] = sr_gmat, er_gmat
         rowblocks_cdiv[0,0], rowblocks_cdiv[1,0] = sr_cdiv, er_cdiv
         rowblocks_dvec[0,0], rowblocks_dvec[1,0] = sr_dvec, er_dvec
         tag,stag = base_tag, base_stag
         for i in np.arange(1,size):
            comm.Recv([rowblocks_gmat[:,i],2,MPI.INT],source=i,tag=tag); tag += 1
            comm.Recv([rowblocks_cdiv[:,i],2,MPI.INT],source=i,tag=tag); tag += 1
            comm.Recv([rowblocks_dvec[:,i],2,MPI.INT],source=i,tag=tag); tag += 1
            
            svec = gmatvals[3*rowblocks_gmat[0,i]:3*(rowblocks_gmat[1,i])]  
            comm.Send([svec,len(svec), MPI.FLOAT], dest=i, tag=stag); stag += 1
            svec = gmatij[:,3*rowblocks_gmat[0,i]:3*(rowblocks_gmat[1,i])].flatten()
            comm.Send([svec,len(svec), MPI.INT], dest=i, tag=stag); stag += 1
            svec = cdiv_temp[rowblocks_cdiv[0,i]:rowblocks_cdiv[1,i]]
            comm.Send([svec,len(svec), MPI.FLOAT], dest=i, tag=stag); stag += 1
            svec = dvec_temp[rowblocks_dvec[0,i]:rowblocks_dvec[1,i]]
            comm.Send([svec,len(svec), MPI.FLOAT], dest=i, tag=stag); stag += 1
         del svec
         gmatvals = gmatvals[3*rowblocks_gmat[0,0]:3*(rowblocks_gmat[1,0])]
         gmatij = gmatij[:,3*rowblocks_gmat[0,0]:3*(rowblocks_gmat[1,0])]
         cdiv_temp = cdiv_temp[rowblocks_cdiv[0,0]:rowblocks_cdiv[1,0]]
         dvec_temp = dvec_temp[rowblocks_dvec[0,0]:rowblocks_dvec[1,0]]

      else:
         tag = base_tag + (rank-1)*3
         comm.Send([np.array([sr_gmat, er_gmat]),2,MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([np.array([sr_cdiv, er_cdiv]),2,MPI.INT], dest=0, tag=tag); tag += 1
         comm.Send([np.array([sr_dvec, er_dvec]),2,MPI.INT], dest=0, tag=tag); tag += 1

         stag = base_stag + (rank-1)*3

         gmatvals = np.empty(3*(er_gmat-sr_gmat), dtype=np.float32)
         gmatij = np.empty(2*3*(er_gmat-sr_gmat), dtype=np.int32)
         cdiv_temp = np.empty((er_gmat-sr_gmat), dtype=np.float32)
         dvec_temp = np.empty((er_gmat-sr_gmat), dtype=np.float32)

         comm.Recv([gmatvals,len(gmatvals),MPI.FLOAT],source=0, tag=stag); stag += 1
         comm.Recv([gmatij,len(gmatij),MPI.INT],source=0, tag=stag); stag += 1
         comm.Recv([cdiv_temp,len(cdiv_temp),MPI.FLOAT],source=0, tag=stag); stag += 1
         comm.Recv([dvec_temp,len(dvec_temp),MPI.FLOAT],source=0, tag=stag); stag += 1
         gmatij = gmatij.reshape(2,-1)

      ### load and assemble gmat, cdiv, and dvec
      for i in np.arange(len(dvec_temp)):
         threei = 3*i
         gmat.setValue(gmatij[0,threei],gmatij[1,threei],gmatvals[threei])
         gmat.setValue(gmatij[0,threei],gmatij[1,threei+1],gmatvals[threei+1])
         gmat.setValue(gmatij[0,threei],gmatij[1,threei+2],gmatvals[threei+2])
         cdiv.setValue(gmatij[0,threei],gmatij[0,threei],cdiv_temp[i])
         dvec.setValue(gmatij[0,threei],dvec_temp[i])
      comm.Barrier()
      dvec.assemblyBegin(); dvec.assemblyEnd()
      gmat.assemblyBegin(assembly=FINAL); gmat.assemblyEnd(assembly=FINAL)
      cdiv.assemblyBegin(assembly=FINAL); cdiv.assemblyEnd(assembly=FINAL)
      comm.Barrier()

      del gmatvals, gmatij, cdiv_temp, dvec_temp

      ### build A and gtilde (=dtilde)
      omega = PETSc.Mat().createAIJ([modlen,p.gen.gmat_rows],comm=comm)
      omega = gmat.matTransposeMult(cdiv,omega)
      cdiv.destroy()

      gtg = PETSc.Mat().createAIJ([modlen,modlen],comm=comm)
      mtilde,gtilde = gtg.getVecs() 
      omega.mult(dvec,gtilde)
      gtg = omega.matMult(gmat,gtg)
      dvec.destroy(); omega.destroy(); gmat.destroy()
      mtilde.set(0)

      if p.gen.cmlambda > 1.e-12:

         PETSc.Sys.Print('building fmat')

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

         fmat = gtg.duplicate(copy=False) 
         fmat.setPreallocationNNZ(nnz=len(fmatval))

         ### get the row addresses for this processor 
         sr_fmat, er_fmat = fmat.getOwnershipRange()

         base_tag, base_stag = 77, 444
         if rank == 0:
            rowblocks_fmat[0,0], rowblocks_fmat[1,0] = sr_fmat, er_fmat
            tag, stag = base_tag, base_stag
            for i in np.arange(1,size):
               comm.Recv([rowblocks_fmat[:,i],2,MPI.INT],source=i,tag=tag); tag += 1
               # get indexes for fmatval and fmatij
               left  = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[0,i]))
               right = np.argmin(np.abs(fmatij[0,:]-rowblocks_fmat[1,i]))
               svec = fmatval[left:right]
               comm.Send([len(svec),1,MPI.INT], dest=i, tag=stag); stag += 1
               comm.Send([svec,len(svec),MPI.FLOAT], dest=i, tag=stag); stag += 1                
               svec = fmatij[:,left:right].flatten()
               comm.Send([svec,len(svec),MPI.INT], dest=i, tag=stag); stag += 1
            fmatval = fmatval[sr_fmat:er_fmat]
            fmatij = fmatij[:,sr_fmat:er_fmat]
         else:
            tag, stag = base_tag + (rank-1), base_stag + (rank-1)*3
            comm.Send([np.array([sr_fmat, er_fmat]),2,MPI.INT], dest=0, tag=tag); tag += 1

            comm.Recv([lenrec,1,MPI.INT], source=0, tag=stag); stag += 1

            fmatval = np.empty(lenrec,dtype=np.float32)
            fmatij  = np.empty(2*lenrec,dtype=np.int32)

            comm.Recv([fmatval,lenrec,MPI.FLOAT], source=0, tag=stag); stag += 1
            comm.Recv([fmatij,2*lenrec,MPI.INT], source=0, tag=stag); stag += 1 
            fmatij = fmatij.reshape(2,-1)

         PETSc.Sys.Print('loading sparse fmat...')
         for i in np.arange(len(fmatval)):
            fmat.setValue(fmatij[0,i],fmatij[1,i],fmatval[i])
         comm.Barrier()
         fmat.assemblyBegin(assembly=FINAL); fmat.assemblyEnd(assembly=FINAL)
         comm.Barrier()
         del fmatval, fmatij

         PETSc.Sys.Print('loading cli...')
         cliv = PETSc.Vec().createMPI(modlen,comm=comm)
         gtg.getDiagonal(cliv)
         cli = gtg.duplicate(copy=False)
         cli.setDiagonal(cliv)
         cliv.destroy()

         PETSc.Sys.Print('building C_m...')
         rhs = gtg.duplicate(copy=False)
         rhs = fmat.matTransposeMult(cli,rhs)
         rhs = rhs.matMult(fmat)
         cli.destroy(); fmat.destroy()

         PETSc.Sys.Print('building A...\n')
         gtg.axpy(p.gen.cmlambda,rhs)
         rhs.destroy()
       

      PETSc.Sys.Print('inverting using PETSc lsqr...') 
      ksp = PETSc.KSP().create(comm)
      ksp.setType('lsqr')
      ksp.pc.setType('none')
      ksp.setOperators(gtg)
      ksp.setFromOptions()
      t0 = time.time()
      ksp.solve(gtilde,mtilde)

      gtg.destroy(); gtilde.destroy(); ksp.destroy()

      PETSc.Sys.Print('lsqr took %10f seconds'%(time.time()-t0))
      PETSc.Sys.Print(' Converged in %d iterations '%(ksp.getIterationNumber()))
      PETSc.Sys.Print(' Tolerance Asked: %e %e %d %d'%(ksp.getTolerances()))
      PETSc.Sys.Print(' Converged Reason: %d'%(ksp.getConvergedReason()))
      PETSc.Sys.Print(' ')

      comm.Barrier() 
      if rank == 0:    
         mtildenump = mtilde[...]; mtilde.destroy()
         self._write_post_model(model=mtildenump,p=p)
      comm.Barrier()

      PETSc._finalize
      MPI.Finalize


 


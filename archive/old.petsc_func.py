
   def PPS(self,p):
      """
      
      Parallel solver implementing a PETSc KSP linear routine over the entire scene as a 
      single patch (PPS = Parallel PETSc Single)

      """

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
 

      print "\nRunning PPS:  Parallel PETSc linear solver on a single patch\n"
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

      print 'building sparse matrices'

      modlen = 3*len(address)
      addind, gmatrow, ott = 0, 0, np.arange(3)
      gmat = PETSc.Mat().createAIJ([p.gen.gmat_rows,modlen],nnz=p.gen.gmat_rows*3,comm=PETSc.COMM_WORLD)
      cdiv = PETSc.Mat().createAIJ([p.gen.gmat_rows,p.gen.gmat_rows],nnz=p.gen.gmat_rows,comm=PETSc.COMM_WORLD)
      dvec = PETSc.Vec().createMPI(p.gen.gmat_rows,comm=PETSc.COMM_WORLD)
      for i in np.arange(p.gen.doline):
         if i % 10 == 0: print i
         if address[addind] >= i*cols or address[addind] < (i+1)*cols: # check if there are good vals in row
            unw,cor,los,cnt = self._read_data_lines(unw,cor,los,dum,fidunw,fidcor,fidlos,p,i)
            for j in np.arange(cols):
               if address[addind] == i*cols+j:
                  aind, jind = 3*addind, 3*j 
                  for m in np.arange(len(p.scenes)):
                     if self._interograte_data(unw=unw[m,j],cor=cor[m,j],los=los[m,jind:jind+3],p=p,m=m):
                        corsq = cor[m,j]**2
                        for k in ott:
                           gmat.setValue(gmatrow, aind+k, los[m,jind+k])
                        dvec.setValue(gmatrow, (unw[m,j]*p.scenes[m].dsp2r))
                        cdiv.setValue(gmatrow, gmatrow, (p.scenes[m].dlooks*corsq/(1.-corsq)))
                        gmatrow += 1
                  addind += 1
                  if addind >= len(address):  addind -= 1
         else:
            self._ffwd_data(fidunw=fidunw,fidcor=fidcor,fidlos=fidlos,p=p,i=i)

      dvec.assemblyBegin(); dvec.assemblyEnd()
      gmat.assemblyBegin(assembly=FINAL); gmat.assemblyEnd(assembly=FINAL)
      cdiv.assemblyBegin(assembly=FINAL); cdiv.assemblyEnd(assembly=FINAL)

      # stow data files and release memory       
      del unw,cor,los,dum
      for m in np.arange(len(p.scenes)):
         fidunw[m].close(); fidcor[m].close(); fidlos[m].close()

      omega = PETSc.Mat().createAIJ([modlen,p.gen.gmat_rows],comm=PETSc.COMM_WORLD)
      omega = gmat.matTransposeMult(cdiv,omega)
      cdiv.destroy()

      gtg = PETSc.Mat().createAIJ([modlen,modlen],comm=PETSc.COMM_WORLD)
      x,b = gtg.getVecs() 
      omega.mult(dvec,b)
      gtg = omega.matMult(gmat,gtg)
      dvec.destroy(); omega.destroy(); gmat.destroy()
      x.set(0)

      if p.gen.cmlambda > 1.e-12:
         print '\nbuilding fmat...'
         fmatt = self._build_laplace_scipy_sparse(a=address,p=p)
         fmat = PETSc.Mat().createAIJ([fmatt.get_shape()[0],fmatt.get_shape()[1]],nnz=fmatt.getnnz(),
                     comm=PETSc.COMM_WORLD)
         fmattinds = fmatt.nonzero()
         print 'loading PETSc fmat...'
         for i in np.arange(fmatt.getnnz()):
            row, col = fmattinds[0][i], fmattinds[1][i]
            fmat.setValue(row, col, fmatt[row,col])
         fmat.assemblyBegin(assembly=FINAL); fmat.assemblyEnd(assembly=FINAL)
         del fmatt, fmattinds

         print 'building cli...'
         cliv = PETSc.Vec().createMPI(modlen,comm=PETSc.COMM_WORLD)
         gtg.getDiagonal(cliv)
         cli = gtg.duplicate(copy=False)
         cli.setDiagonal(cliv)
         cliv.destroy()

         print 'building C_m...'
         rhs = gtg.duplicate(copy=False)
         rhs = fmat.matTransposeMult(cli,rhs)
         rhs = rhs.matMult(fmat)
         cli.destroy(); fmat.destroy()

         print 'building A...'
         gtg.axpy(p.gen.cmlambda,rhs)
         rhs.destroy()
       

      print 'inverting using PETSc...' 
      ksp = PETSc.KSP()
      ksp.create(PETSc.COMM_WORLD)
      ksp.setType('lsqr')
      ksp.getPC().setType('jacobi')
      ksp.setOperators(gtg)
      ksp.setFromOptions()
      t0 = time.time()
      ksp.solve(










 


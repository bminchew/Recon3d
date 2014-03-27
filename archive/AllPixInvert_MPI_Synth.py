import numpy as np
import tsinsar as ts
import sys
import matplotlib.pyplot as plt
import scipy.linalg as lm
import scipy.io as sio
import h5py
import datetime as dt
import mpi4py
from mpi4py import MPI
import petsc4py
petsc4py.init(sys.argv,comm=MPI.COMM_WORLD)
from petsc4py import PETSc

#-------------------------------------------------------
# Communicator
Com = MPI.COMM_WORLD

#-------------------------------------------------------
#-------------------------------------------------------
# Things you should modify

fname = 'DATA/RAW-STACK_truth.h5'			# Input file
oname = 'MP-PYTS-Synth.h5'				# Output file
Conv2h5 = 'no'                				# If you want to store results in h5 file
							# Yes: Stores and sort everything in h5file, This can take a long time, depending on your file system
							# No: Binary output, unsorted, need to use MPfinalize.py

Refx = 75						# x position of ref pixel
Refy = 125						# y position of ref pixel
nref = 10						# Width of the Reference Area in pixels

stax = 0						# First value along x
endx = 350						# Last value along x
stex = 2						# Step along x
stay = 0						# First value along y
endy = 500						# Last value along y
stey = 2						# Step along y

masterind = 0						# Index of the master date

#-------------------------------------------------------
#-------------------------------------------------------

PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print('         Massive Time Series Inversion System          ')
PETSc.Sys.Print(' ') 
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')

#-------------------------------------------------------
# PETSc Stuffs

ZER = PETSc.Mat.Option.IGNORE_ZERO_ENTRIES		# Ignore zero entries in mat() and vec()
INS = PETSc.InsertMode.INSERT_VALUES			# Insert, rather than Add
ADD = PETSc.InsertMode.ADD_VALUES			# Add, rather than insert
Opt = PETSc.Options()

rank = Com.Get_rank()

#-------------------------------------------------------
# Get Informations

PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' Get Design Matrix From the h5file')
PETSc.Sys.Print(' ')

Com.Barrier()

if rank==0:
	Nsar = 90
	fdat = h5py.File(fname,'r')
	fout = h5py.File(oname,'w')
	print 'Process %d: Get Info'%(rank)
	Jmat = fdat['Jmat'].value
	Nifg = Jmat.shape[0]
	Nsar = Jmat.shape[1]

	# Move stuff to the output file
	da = fdat['dates'].value
	ti = fdat['tims'].value
	fout.create_dataset('dates',data=da)
	fout.create_dataset('tims',data=ti)
	fout.create_dataset('masterind',data=masterind)
	
	# Read data
	igram = fdat['figram']
	Ny = igram[:,stay:endy:stey,stax:endx:stex].shape[1]
        Nx = igram[:,stay:endy:stey,stax:endx:stex].shape[2]
	igramsub = fout.create_dataset('Data',(Nifg,Ny,Nx),'f')
        igramsub = igram[:,stay:endy:stey,stax:endx:stex]

	tims = fdat['tims']

	print 'Process %d: Build time function'%(rank)
	rep = [['SBAS',masterind]]
	H,mName,regF = ts.Timefn(rep,tims) 

	print 'Process %d: Build Generic G matrix'%(rank)
	indJ = np.arange(Nsar)
	indJ = np.delete(indJ,masterind)
	indH = np.arange(H.shape[1])
	indH = np.delete(indH,masterind)
	Jmat = Jmat[:,indJ]
	H = H[indJ,:]
	H = H[:,indH]
	nm = H.shape[1]
	fout.create_dataset('Hmat',data=H)

	Gg = np.dot(Jmat,H)
	[xx,yy] = np.where(np.abs(Gg)<1e-20)
	Gg[xx,yy] = 0.0

	print 'Process %d: Build Mask, No tolerance to NaN'%(rank)
	cmask = fdat['cmask']
	nmask = np.sum(igramsub,axis=0)
	cmask = cmask[stay:endy:stey,stax:endx:stex]*nmask
	cmask[np.isfinite(cmask)]=1
	indpixcoh = np.where(np.isfinite(cmask.flatten()))
	yyok,xxok = np.where(np.isfinite(cmask))
	ng = np.size(xxok)
	Npix = np.array(indpixcoh).shape[1]

	Nsar = np.array(Nsar)
	Nifg = np.array(Nifg)
	Ny = np.array(Ny)
	Nx = np.array(Nx)
	nm = np.array(nm)
	Npix = np.array(Npix)
	ng = np.array(ng)
	xxok = np.array(xxok)
	yyok = np.array(yyok)
	cmask = np.array(cmask)
	Gg = np.array(Gg)

else:
	Nsar = None
	Nifg = None
	Ny = None
	Nx = None
	nm = None
	Npix = None
	ng = None
	xxok = None
	yyok = None
	cmask = None
	Gg = None
	Jmat = None

Com.Barrier()

# Broadcast 
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Broadcasting Informations, Mask and Generic G to the %d Processes'%(Com.Get_size()))

Nsar = Com.bcast(Nsar,root=0)
Com.Barrier()
Nifg = Com.bcast(Nifg,root=0)
Com.Barrier()
Ny = Com.bcast(Ny,root=0)
Com.Barrier()
Nx = Com.bcast(Nx,root=0)
Com.Barrier()
nm = Com.bcast(nm,root=0)
Com.Barrier()
Npix = Com.bcast(Npix,root=0)
Com.Barrier()
ng = Com.bcast(ng,root=0)
Com.Barrier()
xxok = Com.bcast(xxok,root=0)
Com.Barrier()
yyok = Com.bcast(yyok,root=0)
Com.Barrier()
cmask = Com.bcast(cmask,root=0)
Com.Barrier()
Gg = Com.bcast(Gg,root=0)
Com.Barrier()
Jmat = Com.bcast(Jmat,root=0)
Com.Barrier()

PETSc.Sys.Print(' Broadcasting is done')
PETSc.Sys.Print(' ')

#-------------------------------------------------------
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Working on %d interferograms from %d SAR acquisitions'%(Nifg,Nsar))
PETSc.Sys.Print(' Size of each image: %d %d'%(Ny,Nx))
PETSc.Sys.Print(' %d data are involved'%(Npix*Nifg))
PETSc.Sys.Print(' ')
#-------------------------------------------------------
#-------------------------------------------------------
# Build and fill a MPI vector for data
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Building a MPI Vector of size %d'%(Npix*Nifg))
PETSc.Sys.Print(' distributed on %d Processe(s)'%(Com.Get_size()))

d = PETSc.Vec().createMPI(Npix*Nifg,comm=Com)
Com.Barrier()

if rank==0:

	# Get Ownerships
	I = []
	I.append(d.getOwnershipRanges())
	I = np.array(I).flatten()

	# Number of Nan per ifg
	NumNan = Nx*Ny - len(indpixcoh)
	
	# Load Sub Igram
	for i in xrange(Nifg):
		Ref = igramsub[i,Refy-(nref/2):Refy+(nref/2),Refx-(nref/2):Refx+(nref/2)]
		Ref = ts.nanmean(Ref.flatten())
		igramsub[i,:,:] = igramsub[i,:,:]*cmask-Ref 
	data = np.swapaxes(np.swapaxes(igramsub,0,1),1,2).flatten()
	ii = np.where(np.isfinite(data))
	data = data[ii]

	# Loop over each Owned Ranges
	for i in xrange(len(I)-1):
		indr = range(I[i],I[i+1])
		if i==0:
			# if you are Proc 0, fill the d vector
			Val = data[indr].flatten()
			d.setValues(indr,Val,INS)
		else:
			# Send data to other Process
			Val = data[indr].flatten()
			#print 'Proc ',rank,' to Proc ',i,' : ',Val
			Com.Send([Val, len(Val), MPI.DOUBLE],dest=i,tag=77)

else:
	# Other Processes receive the data and and fill the d vector	
	Is,Ie = d.getOwnershipRange()
	indr = range(Is,Ie)
	Val = np.empty((Ie-Is,)).astype(np.float64) #,dtype=np.float64)
	Com.Recv([Val, len(Val), MPI.DOUBLE],source=0,tag=77)
	#print 'Proc ',rank,' : ',Val
	d.setValues(indr,Val,INS)

# Wait for number 0 to finish building this guy from data
Com.Barrier()

# Assemble
d.assemblyBegin()
d.assemblyEnd()
# Data vector ready

PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Data Vector ready')
PETSc.Sys.Print(' ')

#-------------------------------------------------------
#-------------------------------------------------------
# Create the solution vector m
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Create the solution vector, m')

# Some size variables
Np = np.int(Npix*nm + 2*(Nsar-1))
Nd = np.int(Npix*Nifg)
nm = np.int(nm)
Nsar = np.int(Nsar)

PETSc.Sys.Print(' Size [%d]'%(Np))
m = PETSc.Vec().createMPI(Np,comm=Com)
m.assemblyBegin()
m.assemblyEnd()

PETSc.Sys.Print(' ')

#-------------------------------------------------------
#-------------------------------------------------------
# Build the G matrix
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')

# Create a sparse G matrix
PETSc.Sys.Print(' Create the G matrix, AIJ Sparse')
PETSc.Sys.Print(' distributed on %d Processe(s)'%(Com.Get_size()))
G = PETSc.Mat().createAIJ([Nd,Np],nnz=[2*nm,2*nm],comm=Com)
gl,gc = G.getSize()
PETSc.Sys.Print(' Size [%d,%d]'%(gl,gc))
PETSc.Sys.Print(' ')
G.setOption(ZER,1)

# Create the Orbit matrix
Orb = np.hstack((Jmat,Jmat))

# Get the OwnershipRange
Istart,Iend = G.getOwnershipRange()

# Check if d Ownership is similar to G Ownership ranges
Isd,Ied = d.getOwnershipRange()
Com.Barrier()
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' While Running, Please check that G and d OwnershipRanges are similar')
PETSc.Sys.Print('           Gstart |  Gend  | dstart |  dend  || Number of Lines')
Com.Barrier()
print 'Proc %3d : %6d | %6d | %6d | %6d || %6d'%(rank,Istart,Iend,Isd,Ied,Iend-Istart)
Com.Barrier()
PETSc.Sys.Print(' ')
PETSc.Sys.Print('-------------------------------------------------------')
Com.Barrier()

# Index for orbits
ico = range(Npix*nm,Np)

t1 = MPI.Wtime()

# To improve the speed of the process that builds G, we start from Istart and go to the next beginning of pixel, then we fill G with big blocks, until there is no room for blocks and we fill the rest
# Istart ---------- Istart+Ifpix --------------------------------- Iend-Ilpix ---------- Iend
#           Line                               Block                             Line
#            by                                  by                               by
#           Line                               Block                             Line 

Ifpix = (Istart/Nifg+1)*Nifg - Istart
Ilpix = Iend - Iend/Nifg*Nifg


PETSc.Sys.Print(' ')
PETSc.Sys.Print('           Istart |  Istart+Ifpix  |  Iend-Ilpix  |  Iend  | Number of Pixels')
Com.Barrier()
print 'Proc %3d : %6d |     %6d     |    %6d    | %6d | %6d '%(rank,Istart,Istart+Ifpix,Iend-Ilpix,Iend,(Iend-Ilpix-Istart-Ifpix)/Nifg)
Com.Barrier()
PETSc.Sys.Print(' ')
PETSc.Sys.Print('-------------------------------------------------------')
Com.Barrier()
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Time Management:')
Com.Barrier()

Com.Barrier()

t11 = MPI.Wtime()
for i in xrange(Istart,Istart+Ifpix):
	
	# Generic Gg
	nbGg = i/Nifg
	RowInGg = i - nbGg*Nifg 
	ic = range(nbGg*nm,nbGg*nm+nm)
	G.setValues(i,ic,Gg[RowInGg,:],INS)

	# Put the Orbits right
	x = (xxok[nbGg].astype(np.float) - np.float(Refy))/(Ny.astype(np.float))
	y = (yyok[nbGg].astype(np.float) - np.float(Refx))/(Nx.astype(np.float))
	O = Orb[RowInGg,:].copy()
	on = np.array(np.where(O!=0.0)).flatten()
	if len(on)==2:
		O[on[0]] = O[on[0]]*x
		O[on[1]] = O[on[1]]*y
	else:
		O[on[0]] = O[on[0]]*x
                O[on[1]] = O[on[1]]*x
		O[on[2]] = O[on[2]]*y
                O[on[3]] = O[on[3]]*y
	G.setValues(i,ico,O,INS)
t12 = MPI.Wtime()
Ft = t12-t11

t11 = MPI.Wtime()
tg = 0
to = 0

for i in xrange(Istart+Ifpix,Iend-Ilpix,Nifg):
	
	# Generic Gg
	nbGg = i/Nifg
	ir = range(nbGg*Nifg,nbGg*Nifg+Nifg)
	ic = range(nbGg*nm,nbGg*nm+nm)
	tt = MPI.Wtime()
	G.setValues(ir,ic,Gg.flatten(),INS)
	ttt = MPI.Wtime()
	tg = tg +(ttt-tt)	

	# Orbits
	x = (xxok[nbGg].astype(np.float) - np.float(Refy))/(Ny.astype(np.float))
        y = (yyok[nbGg].astype(np.float) - np.float(Refx))/(Nx.astype(np.float))
	tt = MPI.Wtime()
	Orb = np.hstack((Jmat*x,Jmat*y))
	G.setValues(ir,ico,Orb.flatten(),INS)
	ttt = MPI.Wtime()
	to = to + (ttt-tt)

t12 = MPI.Wtime()
Mt = t12-t11

# Re-build Orb
Orb = np.hstack((Jmat,Jmat))

t11 = MPI.Wtime()
for i in xrange(Iend-Ilpix,Iend):
	# Generic Gg
        nbGg = i/Nifg
        RowInGg = i - nbGg*Nifg 
        ic = range(nbGg*nm,nbGg*nm+nm)
        G.setValues(i,ic,Gg[RowInGg,:].tolist(),INS)

        # Put the Orbits right
        x = (xxok[nbGg].astype(np.float) - np.float(Refy))/(Ny.astype(np.float))
        y = (yyok[nbGg].astype(np.float) - np.float(Refx))/(Nx.astype(np.float))
        O = Orb[RowInGg,:].copy()
        on = np.array(np.where(O!=0.0)).flatten()
        if len(on)==2:
                O[on[0]] = O[on[0]]*x
                O[on[1]] = O[on[1]]*y
        else:
                O[on[0]] = O[on[0]]*x
                O[on[1]] = O[on[1]]*x
                O[on[2]] = O[on[2]]*y
                O[on[3]] = O[on[3]]*y
        G.setValues(i,ico,O,INS)
t12 = MPI.Wtime()
LT = t12-t11

t2 = MPI.Wtime()
print 'Proc %d | Step 1 : %f | Step 2 : %f %f %f | Step 3 : %f | Total Time : %f '%(rank,Ft,Mt,tg,to,LT,t2-t1)

Com.Barrier()

PETSc.Sys.Print(' Assembling G')
G.assemblyBegin()
G.assemblyEnd()

#-------------------------------------------------------                                                 
#-------------------------------------------------------                                                
# Solve it 
PETSc.Sys.Print('-------------------------------------------------------')                            
PETSc.Sys.Print(' ')                                                                                      

PETSc.Sys.Print(' Initialize KSP Solver')

# Initialize the Solver
ksp = PETSc.KSP().create(comm=Com)

# Default options
ksp.setType('lsqr')
ksp.pc.setType('none')

# Set the Tolerance:
# If you want 0.1 mm tolerance on the whole array, as R=||b-Ax||2
# and the solver stops when R < max(rtol*||b||2,atol)
# you should set atol so that RMS ~= atol/(Npix*Nifg) = 0.1
# Here is an automatic guess, if you give a value in argument, it will erase that guess
#atol = 0.001*np.float(Nd)
#rtol = 1e-6
#ksp.setTolerances(rtol,atol,100000,100000)

# set options from arguments
ksp.setFromOptions()
rtol, atol, it1, it2 = ksp.getTolerances()

## Set the Tolerance:
## set options from arguments
#ksp.setFromOptions()
#ksp.setTolerances(rtol,1e-50,100000,100000)
#rtol, atol, it1, it2 = ksp.getTolerances()

PETSc.Sys.Print(' ')

# set Preconditioner
if ksp.pc.getType()=='none':
	ksp.setOperators(G)
	PETSc.Sys.Print(' Using a %s solver with no preconditioner'%(ksp.getType()))
else:
	PETSc.Sys.Print('********************************************')
	PETSc.Sys.Print('********************************************')
	PETSc.Sys.Print(' ')
	PETSc.Sys.Print('  There is no preconditioner for that case  ')
	PETSc.Sys.Print('  Please write one, that would be cool :-)  ')
	PETSc.Sys.Print(' ')
	PETSc.Sys.Print('********************************************')
	PETSc.Sys.Print('********************************************')
	PETSc._finalize()
	sys.exit(1)

# Print Stuff
PETSc.Sys.Print(' Tolerances asked:  %e %e %d %d'%(rtol, atol, it1, it2))

# Solve!
ksp.solve(d, m)
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Converged in %d iterations '%(ksp.getIterationNumber()))
PETSc.Sys.Print(' Tolerance Asked: %e %e %d %d'%(ksp.getTolerances()))
PETSc.Sys.Print(' Converged Reason: %d'%(ksp.getConvergedReason()))
PETSc.Sys.Print(' ')
ksp.view()
PETSc.Sys.Print(' ')

# Destoy the solver
ksp.destroy()

#-------------------------------------------------------
#-------------------------------------------------------
# Compute the Predicted data vector
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Save to file %s'%(oname))
PETSc.Sys.Print(' Write the Data Residuals')
PETSc.Sys.Print(' ')

# Create Residual Vector
ds = PETSc.Vec().createMPI(Npix*Nifg,comm=Com)

# Multiply G*m
G.mult(m,ds)

# Destroy G
G.destroy()

# Here the writting process is a bit odd, but it works.
# We cannot write to h5 file in parallel because it cannot handle NFS mounted dis, 
# And we want it to work on every machine
# So we write a binary file in parallel, then read it with one process and feed it to the h5 file

# Open the file for output
fres = 'Residu.dat'
fpetscout = PETSc.Viewer()
PETSc.Viewer.createBinary(fpetscout,fres,mode=PETSc.Viewer.Mode.WRITE,format=PETSc.Viewer.Format.ASCII_DENSE,comm=Com)

# Write ds in it
ds.view(viewer=fpetscout)

# Close this guy
fpetscout.destroy()

# Wait until everyone did this
Com.Barrier()

# Destroy these guys
ds.destroy()
d.destroy()

# If you are rank 0
if rank==0:

	if Conv2h5 in ('yes','Yes','y','Y'):

		# Create the output vector in the h5file
		dSim = fout.create_dataset('dSim',(Nifg,Ny,Nx),'f')  
		dSim[:,:,:] = np.nan
	
		# Open the binary file and read the vector
		fpetscin = open(fres,'r')
		ddw = np.fromfile(fpetscin).byteswap()[1:]
		fpetscin.close()
		
		# Create Progress Bar
		print 'Writing for %d pixels'%(len(xxok))
   		toto = ts.ProgressBar(maxValue=len(xxok))
	
		# Loop on xxok and yyok
		for n in range(0,len(xxok)):
			istart = n*Nifg
			iend = istart+Nifg
			dSim[:,yyok[n],xxok[n]] = ddw[istart:iend]
			toto.update(n,every=5)

		# Close Progress Bar
		toto.close()

	else:

		fbinout = open('Xok.dat','w')
		xout = xxok.astype(np.float32)
		xout.tofile(fbinout)
		fbinout.close()
		fbinout = open('Yok.dat','w')
		yout = yyok.astype(np.float32)
		yout.tofile(fbinout)
		fbinout.close()

Com.Barrier()

#-------------------------------------------------------                                                 
#-------------------------------------------------------                                                
# Save the Outputs
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Write the A Posteriori model parameters')
PETSc.Sys.Print(' ')

# Do the same process for m vector
# Open a petsc viewer to store the m vector
fmod = 'Model.dat'
fpetscout = PETSc.Viewer()
PETSc.Viewer.createBinary(fpetscout,fmod,mode=PETSc.Viewer.Mode.WRITE,format=PETSc.Viewer.Format.ASCII_DENSE,comm=Com)

# Write m in Model.dat
m.view(viewer=fpetscout)

# Close the file
fpetscout.destroy()

# Wait for everyone
Com.Barrier()

# Destroy m
m.destroy()

# If you are Rank 0
if rank==0:
	
	if Conv2h5 in ('yes','Yes','y','Y'):

		# Create the input vector in the h5 file
		MPhi = fout.create_dataset('PhiParam',(Nsar-1,Ny,Nx),'f')
		MOrb = fout.create_dataset('OrbParam',(2,nm),'f') 
		MPhi[:,:,:] = np.nan
		MOrb[:,:] = np.nan

		# Open the binary file and read the crap
		fpetscin = open(fmod,'r')
		mw = np.fromfile(fpetscin).byteswap()[1:]

		# Create a progress bar
		toto = ts.ProgressBar(maxValue=len(xxok))
	
		# Loop on non-masked pixels
		for n in range(0,len(xxok)):
			istart = n*(Nsar-1)
			iend = istart+Nsar-1
			MPhi[:,yyok[n],xxok[n]] = mw[istart:iend]
			toto.update(n,every=5)

		# Close the Progress Bar
		toto.close()

		# Store Orbits
		MOrb[0,:] = mw[-2*(Nsar-1):-(Nsar-1)]
		MOrb[1,:] = mw[-(Nsar-1):]

# Wait for everyone
Com.Barrier()

#-------------------------------------------------------                                                 
#-------------------------------------------------------                                                
# Save the Outputs
PETSc.Sys.Print('-------------------------------------------------------')
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Write the Phase Evolution ')
PETSc.Sys.Print(' ')

if rank==0:
	
	if Conv2h5 in ('yes','Yes','y','Y'):

		Phi = fout.create_dataset('Phi',(Nsar,Ny,Nx),'f')
		Orb = fout.create_dataset('Orb',(Nsar,Ny,Nx),'f')
		Phi[:,:,:] = np.nan
		Orb[:,:,:] = np.nan
		Phi[0,:,:] = 0.0
		Orb[0,:,:] = 0.0

		# Open a Progress Bar
		toto = ts.ProgressBar(maxValue=Nx*Ny)
		uc = 0
	
		for j in range(Ny):
			for i in range(Nx):
				Phi[1:,j,i] = np.dot(H,MPhi[:,j,i])
				ii = (np.float(i) - np.float(Refx))/np.float(Nx)
				jj = (np.float(j) - np.float(Refy))/np.float(Ny)
				o = np.ones((MOrb.shape))
				o[0,:] = ii
				o[1,:] = jj
				Orb[1:,j,i] = MOrb[0,:]*o[0,:]+MOrb[1,:]*o[1,:]
				uc = uc + 1
				toto.update(uc,every=5)

		# Close the Progress Bar
		toto.close()

	else:

		print 'To sort the data in the newly created dat files'
		print '   Copy your h5 output file to a fast i/o disk space, run MPIfinalize.py and get back the h5file:'
		print '   rm -rf /scratch/jolivetr/%s'%oname
		print '   cp %s /scratch/jolivetr/.'%oname
		print '   MPfinalize.py -i %d -s %d -x %d -y %d -refx %d -refy %d -o /scratch/jolivetr/%s '%(Nifg,Nsar,Nx,Ny,Refx,Refy,oname)
		print '   cp /scratch/jolivetr/%s . '%oname
		print ' '

#-------------------------------------------------------
#-------------------------------------------------------
# Finalize the Job
PETSc.Sys.Print('-------------------------------------------------------') 
PETSc.Sys.Print(' ')
PETSc.Sys.Print(' Finalize the Job')
PETSc.Sys.Print(' ')
PETSc.Sys.Print('-------------------------------------------------------') 
# Close h5file
if rank==0:
	fdat.close()
	fout.close()

# Finalize PETSc
PETSc._finalize

# Finalize MPI
MPI.Finalize

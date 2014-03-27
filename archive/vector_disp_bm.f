   program vector_disp

!****************************************************************



!****************************************************************
!**     
!**   FILE NAME: vector_disp_bm.f
!**     
!**   DATE WRITTEN: September 2011
!**     
!**   PROGRAMMER: 
!**                          Brent Minchew  
!**                California Institute of Technology 
!**                       bminchew@caltech.edu
!**      
!**   FUNCTIONAL DESCRIPTION: 3D vector inversion estimated using
!**               weighted least squares from unwrapped phase or
!**               line of sight displacement files
!**
!**   INPUTS:  Phase or LOS displacement, correlation, LOS file
!**            (LOS must be pixel-by-pixel interleave in order ENU)
!**
!**   OUTPUTS: 1) Velocity field in component order and reference
!**                  frame as given in LOS file
!**            2) [optional] error estimation for each component
!**            3) [optional] mean square error for each pixel
!**   
!**   ROUTINES CALLED:  dggglm, dgetri, dgetrf
!**     
!**   NOTES: must be called with appropriate command file 
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed                  
!**   ------------       ----------------             
!** Jan 8, 2013   Replaced MIGS inversion with LAPACK routines: dggglm, dgetri
!**     
!*****************************************************************


   implicit none

   character*250           cmdfile,otpref,otfold,eastnm,nornm,upnm,msenm,gtgnm
   character*250           str,str2,label,dumname,vmagnm,numnm,acornm,gdopnm
   character*1             errstr 
   character(len=250),dimension(500)::  unwnm,cornm,losnm,aznm

   integer                 ii,jj,kk,mm,zz,ierr,ios,ar,r,k,ns,iid,dumi,tp,tf,alaz
   integer                 numscenes,master,cnt,doline,errest,nonul,maxcol,vm
   integer                 stv(13),today(3),now(3),indx(3),mseout,scncnt,eoc,getwdp
   integer                 info,lwork,nrhs,corout,gdopest,gtgout,gtgoff,ocoff
   real*4                  lambda,timeconv,distconv,pi,dum1,dum2,lne,nulo
   integer,allocatable::   lines(:),cols(:),colsrt(:),azlogic(:),ipiv(:)
   real*4,allocatable::    unw(:,:),cor(:,:),los(:,:),east(:),north(:),up(:)
   real*4,allocatable::    dumv1(:),dumv2(:),azoff(:,:),ehdg(:),nhdg(:),hdg(:)
   real*4,allocatable::    erre(:),errn(:),erru(:),nulaz(:),gdopvec(:)
   real*4,allocatable::    latin(:),lonin(:),latspc(:),lonspc(:),nulu(:),nulc(:),nullos(:)
   real*4,allocatable::    nlooks(:),r2dsp(:),msev(:),vmag(:),numout(:),avecor(:)
   real*4,allocatable::    gtgeast(:),gtgnorth(:),gtgup(:),gtgoen(:),gtgonu(:),gtgoeu(:)
   real*4,allocatable::    ocoen(:),oconu(:),ocoeu(:),wdop(:)
   real*8,allocatable::    dgmat(:,:),ddvec(:),dvmat(:,:),dbmat(:,:),dxvec(:),dyvec(:)
   real*8,allocatable::    ccoef(:),work(:)
   real*4                  gtg(3,3),gtgc(3,3)
   real*8                  dgtg(3,3),gdopmat(3,3),mse,dcor,dscor
   integer,allocatable::   scol(:),acol(:),bcol(:),latoff(:),lonoff(:),eline(:),sline(:)
   integer,allocatable::   elinesrt(:),logic(:)
   integer,allocatable::   losid(:),unwid(:),corid(:),azid(:)     ! file ids 

   type :: parstr
      real*4      pi,fpii,r2d
   end type
   type(parstr) par

   character(len=11),parameter:: fu='unformatted'
   character(len=6), parameter:: di='direct'
   character(len=7), parameter:: re='replace'


   print*,' ';print*,' '
   if (iargc().eq.1) goto 900
      print*,'  ~~~  3D Displacement Inversion  ~~~  '
      print*,' '
      print*,' Usage vector_disp_bm input_command_file '
      print*,' '
      print*,' '
      stop

900  continue

   call getarg(1,cmdfile)
   open(3,file=cmdfile,status='old')

   ios         = 0
   lne         = 1.d0    
   errest      = 0
   timeconv    = 1.d0
   distconv    = 1.d0
   pi          = 4.d0*atan(1.d0)
   numscenes   = 5000
   ns          = 0
   mseout      = 0
   eoc         = 0
   vm          = 0
   alaz        = 0
   corout      = 0
   gdopest     = 0
   gtgout      = 0
   gtgoff      = 0
   ocoff       = 0
   getwdp      = 0

   scncnt = 0

   do while (ios.eq.0.and.ns.le.numscenes.and.scncnt.le.numscenes.and.eoc.eq.0)
      read(3,'(A)',iostat=ios) str
      k = scan(str,'::')
      label = trim(adjustl(str(1:k-1)))
      str=str(k+2:)

      select case (label)

      case('Number of Scenes (-)')
         read(str,*,iostat=ios) numscenes
         if (numscenes.lt.3) then
            print*,'You must have at least 3 scenes to do the inversion'
            stop
         endif
         allocate(lines(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lines'
         allocate(cols(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in cols'
         allocate(latin(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latin'
         allocate(lonin(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonin'
         allocate(latspc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latspc'
         allocate(lonspc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonspc'
         allocate(nlooks(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nlooks'
         allocate(nulc(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulc'
         allocate(nulu(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulu'
         allocate(nullos(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nullos'
         allocate(nulaz(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nulaz'
         allocate(r2dsp(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in r2dsp'
         allocate(azlogic(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in azlogic'
         allocate(ehdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in ehdg'
         allocate(nhdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in nhdg'
         allocate(hdg(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in hdg'
         azlogic     = 0

      case('Master Scene (-)')
         read(str,*,iostat=ios) master
      case('Radar Wavelength (in preferred output distance units)')
         read(str,*,iostat=ios) lambda
      case('Time converstion factor (applied to all scenes equally)')
         if (trim(adjustl(str)).ne.'none') read(str,*,iostat=ios) timeconv
      case('Distance converstion factor (applied to all scenes equally)')
         if (trim(adjustl(str)).ne.'none') read(str,*,iostat=ios) distconv

      case('Output file prefix')
         otpref = str
      case('Output folder')
         otfold = str
      case('Output error estimates')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output diag(G^T*W*G)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output model covariance (G^TWG)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            errest = 1
         endif
      case('Output offdiag(G^T*W*G)^-1)')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            ocoff = 1
         endif
      case('Output diag(G^TG)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgout = 1
         endif
      case('Output diag(G^T*G)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgout = 1
         endif
      case('Output offdiag(G^T*G)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgoff = 1
         endif
      case('Output offdiag(G^TG)^-1')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gtgoff = 1
         endif
      case('Output sqrt(tr((G^T*W*G)^-1))')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            getwdp = 1
         endif
      case('Output GDOP sqrt(tr((G^TG)^-1))')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            gdopest = 1
         endif
      case('Output mean square error estimate')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            mseout = 1
         endif
      case('Output velocity magnitude')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            vm = 1
         endif
      case('Output average correlation')
         errstr = trim(adjustl(str))
         if (errstr.eq.'y'.or.errstr.eq.'Y') then
            corout = 1
         endif

      case('Output null value')
         if (trim(adjustl(str)).eq.'NaN'.or.trim(adjustl(str)).eq.'nan') then
            lne  = 0.d0
            nulo = 0.d0
         elseif (trim(adjustl(str)).eq.'Inf'.or.trim(adjustl(str)).eq.'inf') then
            lne  = 0.d0
            nulo = 1.d0
         elseif (trim(adjustl(str)).eq.'-Inf'.or.trim(adjustl(str)).eq.'-inf') then
            lne  = 0.d0
            nulo = -1.d0
         else
            read(str,*,iostat=ios) nulo
            lne = 1.d0
         endif
         nulo = nulo/lne

      case('Scene')
         read(str,*,iostat=ios) ns
      case('Samples in Unwrapped or LOS displacement (-)')
         if (ns.gt.0) then
            read(str,*,iostat=ios) cols(ns) 
            scncnt = scncnt + 1
         endif
      case('Unwrapped or LOS displacement file')
         if (ns.gt.0) unwnm(ns) = str
      case('Type of unwrapped or LOS displacement file (enter phase or disp)')
         if (trim(adjustl(str)).eq.'phase'.and.ns.le.numscenes.and.ns.gt.0) then
            r2dsp(ns) = lambda/(4.d0*pi)*timeconv*distconv
         elseif (trim(adjustl(str)).eq.'disp'.and.ns.le.numscenes.and.ns.gt.0) then
            r2dsp(ns) = timeconv*distconv
         elseif (ns.gt.numscenes.or.ns.le.0) then
            continue
         else
            print *,'Invalid type of displacement file in scene ',ns
            print *,'Stopping...'
            stop
         endif
      case('Correlation file')
         if (ns.gt.0) cornm(ns) = str
      case('LOS file (must be by-pixel interleave in order ENU)')
         if (ns.gt.0) losnm(ns) = str
      case('Azimuth offset file (enter none if unavailable)')
         if  (trim(adjustl(str)).ne.'none'.and.trim(adjustl(str)).ne.'None' &
               & .and.trim(adjustl(str)).ne.'NONE'.and.ns.gt.0) then
            aznm(ns) = str
            azlogic(ns) = 1   
            alaz = 1
         endif
      case('Peg heading (only needed if azimuth offset file is provided)')
         if(azlogic(ns).eq.1) read(str,*,iostat=ios) hdg(ns)
      case('Azimuth offset null value')
         if(azlogic(ns).eq.1) read(str,*,iostat=ios) nulaz(ns)
      case('Upper left corner Latitude (deg)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) latin(ns)
      case('Upper left corner Longitude (deg)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) lonin(ns)
      case('Latitude Spacing (deg/pixel)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) latspc(ns)
      case('Longitude Spacing (deg/pixel)')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) lonspc(ns)
      case('Correlation null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulc(ns)
      case('Phase or displacement null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nulu(ns)
      case('LOS null value')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nullos(ns)
      case('Number of looks')
         if(ns.le.numscenes.and.ns.gt.0) read(str,*,iostat=ios) nlooks(ns)
         
      
      case('End of command file')
         eoc = 1

      case default
         continue
      end select

   enddo

   allocate(colsrt(numscenes))
   colsrt = cols
   call qsort(colsrt,numscenes) 
   maxcol = colsrt(numscenes)
   deallocate(colsrt)
!  now allocate the rest
   allocate(acol(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in acol'
   allocate(bcol(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in bcol'
   allocate(scol(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in scol'
   allocate(eline(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in eline'
   allocate(sline(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in sline'
   allocate(latoff(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in latoff'
   allocate(lonoff(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in lonoff'
   allocate(losid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in losid'
   allocate(unwid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unwid'
   allocate(corid(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in corid'
   allocate(ccoef(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in ccoef'
   allocate(elinesrt(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in elinesrt'
   allocate(logic(numscenes),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in logic'
   allocate(dumv1(10*maxcol),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in dumv1'
   allocate(dumv2(maxcol),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in dumv2'

   allocate(unw(numscenes,cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in unw'
   allocate(cor(numscenes,cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in cor'
   allocate(los(numscenes,3*cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in los'
   allocate(east(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in east'
   allocate(north(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in north'
   allocate(up(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in up'
   allocate(erre(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in erre'
   allocate(errn(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in errn'
   allocate(erru(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in erru'
   allocate(numout(cols(master)),stat=ierr)
            if(ierr.ne.0) print*,'allocation error in numout'
   if (mseout.eq.1) allocate(msev(cols(master)))
   if (vm.eq.1) allocate(vmag(cols(master)))    
   if (gdopest.eq.1) allocate(gdopvec(cols(master)))        
   if (corout.eq.1) allocate(avecor(cols(master)))
   if (getwdp.eq.1) allocate(wdop(cols(master)))
   if (gtgout.eq.1) then
      allocate(gtgeast(cols(master)))
      allocate(gtgnorth(cols(master)))
      allocate(gtgup(cols(master)))
   endif
   if (ocoff.eq.1) then
      allocate(ocoen(cols(master)))
      allocate(oconu(cols(master)))
      allocate(ocoeu(cols(master)))
   endif
   if (gtgoff.eq.1) then
      allocate(gtgoen(cols(master)))
      allocate(gtgonu(cols(master)))
      allocate(gtgoeu(cols(master)))
   endif
   if (alaz.eq.1) then
      allocate(azoff(numscenes,cols(master)))
      allocate(azid(numscenes))
   endif

!  open output files
   otpref = adjustl(otpref)
   otfold = adjustl(otfold)
   tp = len_trim(otpref)
      if (otpref(tp:tp).eq.'.') tp = tp-1
   tf = len_trim(otfold)
      if (otfold(tf:tf).eq.'/') tf = tf-1
   eastnm = otfold(1:tf)//'/'//otpref(1:tp)//'.east'
   nornm  = otfold(1:tf)//'/'//otpref(1:tp)//'.north'
   upnm   = otfold(1:tf)//'/'//otpref(1:tp)//'.up'
   numnm  = otfold(1:tf)//'/'//otpref(1:tp)//'.num'

   open(11,file=eastnm,access=di,form=fu,status=re,recl=4*cols(master))
   open(12,file=nornm,access=di,form=fu,status=re,recl=4*cols(master))
   open(13,file=upnm,access=di,form=fu,status=re,recl=4*cols(master))
   open(19,file=numnm,access=di,form=fu,status=re,recl=4*cols(master))

   if (errest.eq.1) then
      eastnm = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.east'
      nornm  = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.north'
      upnm   = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.up'
      open(14,file=eastnm,access=di,form=fu,status=re,recl=4*cols(master))
      open(15,file=nornm,access=di,form=fu,status=re,recl=4*cols(master))
      open(16,file=upnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif
   if (mseout.eq.1) then
      msenm = otfold(1:tf)//'/'//otpref(1:tp)//'.mse'
      open(17,file=msenm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (vm.eq.1) then
      vmagnm = otfold(1:tf)//'/'//otpref(1:tp)//'.mag'
      open(18,file=vmagnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (corout.eq.1) then
      acornm = otfold(1:tf)//'/'//otpref(1:tp)//'.avecor'
      open(20,file=acornm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (gdopest.eq.1) then
      gdopnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gdop'
      open(21,file=gdopnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (gtgout.eq.1) then
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.east'
      open(22,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.north'
      open(23,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.up'
      open(24,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (gtgoff.eq.1) then
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.en'
      open(25,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.nu'
      open(26,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.gtgi.eu'
      open(27,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (ocoff.eq.1) then
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.en'
      open(28,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.nu'
      open(29,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.eu'
      open(30,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif

   if (getwdp.eq.1) then
      gtgnm = otfold(1:tf)//'/'//otpref(1:tp)//'.obscov.dop'
      open(31,file=gtgnm,access=di,form=fu,status=re,recl=4*cols(master))
   endif 

!  set up each scene
   iid = 41     ! starting input file id 

   do ii=1,numscenes

       unwid(ii) = iid
931      iid = iid + 1
         if (iid.ge.900.and.iid.le.940) goto 931
         if (iid.ge.100.and.iid.le.110) goto 931
       corid(ii) = iid
932      iid = iid + 1
         if (iid.ge.900.and.iid.le.940) goto 932
         if (iid.ge.100.and.iid.le.110) goto 932
       losid(ii) = iid
933      iid = iid + 1
         if (iid.ge.900.and.iid.le.940) goto 933
         if (iid.ge.100.and.iid.le.110) goto 933
       

      dumname = trim(adjustl(unwnm(ii)))
       open(unwid(ii),file=dumname,access=di,form=fu,status='old',recl=4*cols(ii))
       call stat(dumname,stv)
         lines(ii) = stv(8)/(4*cols(ii))
      dumname = trim(adjustl(cornm(ii)))
       open(corid(ii),file=dumname,access=di,form=fu,status='old',recl=4*cols(ii))
      dumname = trim(adjustl(losnm(ii)))
       open(losid(ii),file=dumname,access=di,form=fu,status='old',recl=3*4*cols(ii))
      if (azlogic(ii).eq.1) then
         azid(ii) = iid
934      iid = iid + 1
         if (iid.ge.900.and.iid.le.940) goto 934
         if (iid.ge.100.and.iid.le.110) goto 934
         
         dumname = trim(adjustl(aznm(ii)))
         open(azid(ii),file=dumname,access=di,form=fu,status='old',recl=4*cols(ii))
         azoff(ii,:) = nulaz(ii)

         ehdg(ii) = sin(hdg(ii)*pi/180.d0)
         nhdg(ii) = cos(hdg(ii)*pi/180.d0)
      endif
         

!  calculate pixel offsets for each scene relative to the master
            dumi = 0
       dum1 = (latin(master) - latin(ii))/latspc(ii)
            if (abs(dum1-int(dum1)) .ge. 0.5d0) dumi = dum1/abs(dum1)       
       latoff(ii) = -(int(dum1) + dumi)     ! positive southward
            dumi = 0
       dum1 = (lonin(master) - lonin(ii))/lonspc(ii)
            if (abs(dum1-int(dum1)) .ge. 0.5d0) dumi = dum1/abs(dum1)
       lonoff(ii) = -(int(dum1) + dumi)     !  positive eastward

       eline(ii) = lines(ii) + latoff(ii)
       sline(ii) = -latoff(ii) 
   
!  set column offset conditions 
      if (lonoff(ii).le.0) then      ! shifted left or not at all relative to master
         bcol(ii) = 0                ! no buffer on left
         scol(ii) = -lonoff(ii) + 1  ! starting column
         dum1 = cols(ii) - scol(ii) + 1
            if (dum1.lt.cols(master)) then   !  needs null-valued buffer on the right
               acol(ii) = cols(master) - dum1
            else                             !  no buffer on the right
               acol(ii) = 0
            endif
      else                       !  shifted right relative to master
         bcol(ii) = lonoff(ii)   ! number of null-valued columns to buffer on left
         scol(ii) = 1
         dum1 = cols(ii)
         if (dum1.lt.cols(master)) then
            acol(ii) = cols(master) - dum1   ! number of null-valued columns to buffer onthe right
         else
            acol(ii) = 0                     ! no buffer on the right
         endif
      endif
      ccoef(ii) = 2.d0*nlooks(ii)
   enddo

   elinesrt = eline
   call qsort(elinesrt,numscenes)

   if (lines(master).le.elinesrt(numscenes-2)) then
      doline = lines(master)
   else
      doline = elinesrt(numscenes-2)
   endif

   print*,' '
   print*,' ~~~   Weighted least squares velocity inversion   ~~~'
   print*,' '
   print*,' Total lines    = ',doline
   print*,' Total samples  = ',cols(master)
   print*,' Total scenes   = ',numscenes 
   print*,' ';print*,' '

!  Let's do it
   dumv2 = 0.d0 
   do ii=1,doline

      if(mod(ii,250).eq.0.or.ii.eq.doline.or.ii.eq.1) then
         write(*,"(A)",advance='no') '\b\b\b\b\b\b\b\b\b\b\b\b\b\b'
         write(*,100,advance='no') 'Line = ',ii
      endif

      cnt = 0
      do mm=1,numscenes

         if (ii.ge.latoff(mm).and.ii.le.eline(mm).and.(sline(mm)+ii).gt.0) then
            dumv1 = nulu(mm)
            read(unwid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
               unw(mm,:) = dumv1(scol(mm):cols(master)+scol(mm)-1)
            dumv1 = nulc(mm)
            read(corid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm)) 
               cor(mm,:) = dumv1(scol(mm):cols(master)+scol(mm)-1)
            dumv1 = nullos(mm)
            read(losid(mm),rec=ii+sline(mm)) dumv1(3*bcol(mm)+1:3*(bcol(mm)+cols(mm))) 
               los(mm,:) = dumv1(3*(scol(mm)-1)+1:3*(cols(master)+scol(mm)-1))
            if(azlogic(mm).eq.1) then
               dumv1 = nulaz(mm)
               read(azid(mm),rec=ii+sline(mm)) dumv1(bcol(mm)+1:bcol(mm)+cols(mm))
               azoff(mm,:) = dumv1(scol(mm):cols(master)+scol(mm)-1)
            endif

            cnt = cnt + 1
         endif
      enddo

      east  = nulo
      north = nulo
      up    = nulo 

      erre  = nulo      
      errn  = nulo
      erru  = nulo

      numout = nulo

      if (mseout.eq.1) msev = nulo
      if (vm.eq.1) vmag = nulo
      if (corout.eq.1) avecor = nulo
      if (gdopest.eq.1) gdopvec = nulo
      if (getwdp.eq.1) wdop = nulo
      if (gtgout.eq.1) then
         gtgeast  = nulo
         gtgnorth = nulo
         gtgup    = nulo
      endif
      if (gtgoff.eq.1) then
         gtgoen = nulo
         gtgonu = nulo
         gtgoeu = nulo
      endif    
      if (ocoff.eq.1) then
         ocoen = nulo
         oconu = nulo
         ocoeu = nulo
      endif      

      if (cnt.ge.3) then
         do jj=1,cols(master)
            nonul = 0        
            logic = 0

            do mm=1,numscenes
               dum1 = product(los(mm,(1+(jj-1)*3):(3+(jj-1)*3))-nullos)*cor(mm,jj) ! test null vals
               if(unw(mm,jj).ne.nulu(mm).and.dum1.ne.0.d0.and.cor(mm,jj).ne.nulc(mm)) then
                  nonul = nonul + 1
                  logic(mm) = 1
                  if (azlogic(mm).eq.1.and.azoff(mm,jj).ne.nulaz(mm)) nonul = nonul + 1
               endif
            enddo
      
            numout(jj) = nonul

            if (nonul.ge.3) then
               lwork = 3+nonul+nonul
               allocate(dgmat(nonul,3))
               allocate(dvmat(nonul,nonul))
               allocate(dbmat(nonul,nonul))
               allocate(ddvec(nonul))
               allocate(dxvec(3))
               allocate(dyvec(nonul))
               allocate(work(lwork))

               dvmat = 0.d0
               dbmat = 0.d0 
               dscor = 0.d0
               kk = 1
               do mm=1,numscenes
                  if (logic(mm).eq.1) then
                     dcor = (dble(cor(mm,jj)))**2.d0
                     dvmat(kk,kk)= ccoef(mm)*dcor/(1-dcor)
                     dbmat(kk,kk)= 1.d0/dvmat(kk,kk)
                     dgmat(kk,1) = dble(los(mm,(1+(jj-1)*3)))
                     dgmat(kk,2) = dble(los(mm,(2+(jj-1)*3)))
                     dgmat(kk,3) = dble(los(mm,(3+(jj-1)*3)))
                     ddvec(kk)   = dble(unw(mm,jj)*r2dsp(mm))
                     dscor = dscor + dble(cor(mm,jj))

                     kk = kk + 1
                     if(azlogic(mm).eq.1.and.azoff(mm,jj).ne.nulaz(mm)) then
                        dvmat(kk,kk) = ccoef(mm)*dcor/(1-dcor)
                        dbmat(kk,kk) = 1.d0/dvmat(kk,kk)
                        dgmat(kk,1)  = dble(ehdg(mm))
                        dgmat(kk,2)  = dble(nhdg(mm))
                        dgmat(kk,3)  = 0.d0
                        ddvec(kk)    = dble(azoff(mm,jj))
                        dscor = dscor + dble(cor(mm,jj))
                        kk = kk + 1
                     endif
                  endif
               enddo
               
               dgtg = matmul(matmul(transpose(dgmat),dvmat),dgmat)
               gdopmat = matmul(transpose(dgmat),dgmat)

               call dggglm(nonul,3,nonul,dgmat,nonul,dbmat,nonul,ddvec,dxvec,dyvec,work,lwork,info)
                  ! dggglm is the LAPACK general least squares routine
               east(jj)  = sngl(dxvec(1))
               north(jj) = sngl(dxvec(2))
               up(jj)    = sngl(dxvec(3))

               if (corout.eq.1) avecor(jj) = sngl(dscor/dble(nonul))
               if (vm.eq.1) vmag(jj) = sngl(dsqrt(dxvec(1)**2+dxvec(2)**2+dxvec(3)**2))
               if (gtgout.eq.1.or.gdopest.eq.1.or.gtgoff.eq.1) then
                  call dinv(3,gdopmat)
                  if (gtgout.eq.1) then
                     gtgeast(jj)  = sngl(gdopmat(1,1))
                     gtgnorth(jj) = sngl(gdopmat(2,2))
                     gtgup(jj)    = sngl(gdopmat(3,3))
                  endif
                  if (gtgoff.eq.1) then
                     gtgoen(jj) = sngl(gdopmat(1,2))
                     gtgonu(jj) = sngl(gdopmat(2,3))
                     gtgoeu(jj) = sngl(gdopmat(1,3))
                  endif
                  if (gdopest.eq.1) gdopvec(jj) = sngl(dsqrt(gdopmat(1,1) + gdopmat(2,2) + gdopmat(3,3)))
               endif
               if (errest.eq.1.or.ocoff.eq.1.or.getwdp.eq.1) then
                  call dinv(3,dgtg)  
                  if (errest.eq.1) then
                     erre(jj) = sngl(dgtg(1,1))
                     errn(jj) = sngl(dgtg(2,2))
                     erru(jj) = sngl(dgtg(3,3))
                  endif
                  if (ocoff.eq.1) then
                     ocoen(jj) = sngl(dgtg(1,2))
                     oconu(jj) = sngl(dgtg(2,3))
                     ocoeu(jj) = sngl(dgtg(1,3))
                  endif
                  if (getwdp.eq.1) wdop(jj) = sngl(dsqrt(dgtg(1,1) + dgtg(2,2) + dgtg(3,3)))
               endif
               if (mseout.eq.1) then 
                     if (nonul.gt.2) then
                        mse = sum(dyvec**2)/dble(nonul)
                     else
                        mse = 0.d0
                     endif
                     msev(jj) = sngl(mse)
               endif

               deallocate(dgmat,dvmat,dbmat,ddvec,work,dxvec,dyvec) 
            endif
         enddo
      endif 

      write(11,rec=ii) (east(mm),mm=1,cols(master))
      write(12,rec=ii) (north(mm),mm=1,cols(master))
      write(13,rec=ii) (up(mm),mm=1,cols(master))
      write(19,rec=ii) (numout(mm),mm=1,cols(master))

      if (errest.eq.1) then
         write(14,rec=ii) (erre(mm),mm=1,cols(master))
         write(15,rec=ii) (errn(mm),mm=1,cols(master))
         write(16,rec=ii) (erru(mm),mm=1,cols(master))
      endif
      if (mseout.eq.1) write(17,rec=ii) (msev(mm),mm=1,cols(master))
      if (vm.eq.1) write(18,rec=ii) (vmag(mm),mm=1,cols(master))
      if (corout.eq.1) write(20,rec=ii) (avecor(mm),mm=1,cols(master))
      if (gdopest.eq.1) write(21,rec=ii) (gdopvec(mm),mm=1,cols(master))
      if (gtgout.eq.1) then
         write(22,rec=ii) (gtgeast(mm),mm=1,cols(master))
         write(23,rec=ii) (gtgnorth(mm),mm=1,cols(master))
         write(24,rec=ii) (gtgup(mm),mm=1,cols(master))
      endif
      if (gtgoff.eq.1) then
         write(25,rec=ii) (gtgoen(mm),mm=1,cols(master))
         write(26,rec=ii) (gtgonu(mm),mm=1,cols(master))
         write(27,rec=ii) (gtgoeu(mm),mm=1,cols(master))
      endif 
      if (ocoff.eq.1) then
         write(28,rec=ii) (ocoen(mm),mm=1,cols(master))
         write(29,rec=ii) (oconu(mm),mm=1,cols(master))
         write(30,rec=ii) (ocoeu(mm),mm=1,cols(master))
      endif
      if (getwdp.eq.1) write(31,rec=ii) (wdop(mm),mm=1,cols(master))

   enddo

   print*,' ';print*,' '
   print*,'Done!'
   print*,' ';print*,' '


100   format(a8,i6)
110   format(a15,f9.2)

   end program


!  *************************************************************************************
!  *************************************************************************************

   subroutine dinv(n,mat)
!  returns inv(mat) for a double precision nXn matrix mat 
   implicit none
   integer     n,info,ipiv(n)
   real*8      mat(n,n), work(n*n)
   call dgetrf(n,n,mat,n,ipiv,info)
   call dgetri(n,mat,n,ipiv,work,n*n,info)
   end subroutine


!  *************************************************************************************
!  *************************************************************************************

   subroutine getmse(v,b,g,gtg,m,n,mse)
!  calculates mean square error scalar mse
!
!  inputs:   v = inverse covariance matrix 
!            b = data vector
!            g = linear operator 
!            gtg = inv(transpose(g)*inv(v)*g)
!            m = number of observations
!            n = size of model (will be 3 for 3D vector field)

   implicit none

   integer  m,n
   real*4   gtg(n,n),v(m,m),g(m,n),b(m)
   real*4   va(m,n),vagtg(m,n),at(n,m),atv(n,m)
   real*4   big(m,m),par1(m),mse

   if (m.eq.n) then
      mse = 0
      return
   endif

   at = transpose(g)
   va = matmul(v,g)
   vagtg = matmul(va,gtg)
   atv = matmul(at,v)
   big = v - matmul(vagtg,atv)

   par1(1) = b(1)*big(1,1) + b(2)*big(2,1) + b(3)*big(3,1)
   par1(2) = b(1)*big(1,2) + b(2)*big(2,2) + b(3)*big(3,2)
   par1(3) = b(1)*big(1,3) + b(2)*big(2,3) + b(3)*big(3,3)
 
   mse  = (par1(1)*b(1) + par1(2)*b(2) + par1(3)*b(3))/(m-n)

   end subroutine
   
      


!  *************************************************************************************
!  *************************************************************************************


      SUBROUTINE qsort(COUNT,N)

!     Quicksort   
!C     .. Scalar Arguments ..
      INTEGER N
!C     ..
!C     .. Array Arguments ..
      REAL COUNT(*)
!C     ..
!C     .. Local Scalars ..
      REAL AV,X
      INTEGER I,IF,IFEND,IFK,IFKA,IK,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,MLOOP
!C     ..
!C     .. Local Arrays ..
      INTEGER MARK(100)
!C     ..
!C     .. Executable Statements ..
!C  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED
      IF (N.EQ.1) GO TO 280
      IF (N.GE.1) GO TO 110
      WRITE (6,FMT=100) 

100 FORMAT (/,/,/,20X,' ***KB05A*** ','NO NUMBERS TO BE SORTED ** RETURN TO CALLING PROGRAM')

      GO TO 280
!C  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
!C  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
  110 M = 12
!C  SET UP INITIAL VALUES.
      LA = 2
      IS = 1
      IF = N
      DO 270 MLOOP = 1,N
!C  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 140
!C********* FINAL SORTING ***
!C  ( A SIMPLE BUBBLE SORT )
        IS1 = IS + 1
        DO 130 J = IS1,IF
          I = J
  120     IF (COUNT(I-1).LE.COUNT(I)) GO TO 130
          AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          I = I - 1
          IF (I.GT.IS) GO TO 120
  130   CONTINUE
        LA = LA - 2
        GO TO 260 
!C             *******  QUICKSORT  ********
!C  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
!C  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
!C  HIGHEST ADDRESS.
  140   IY = (IS+IF)/2
        X = COUNT(IY)
        COUNT(IY) = COUNT(IF)
!C  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
!C  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
!C  OF X .
        K = 1
        IFK = IF
!C  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
!C  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AS NECESSARY,
!C  UNTIL THEY MEET .
        DO 160 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 160
          IF (I.GE.IFK) GO TO 170
          COUNT(IFK) = COUNT(I)
          K1 = K
          DO 150 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GE.X) GO TO 150
            IF (I.GE.IFK) GO TO 180
            COUNT(I) = COUNT(IFK)
            GO TO 160

  150     CONTINUE
          GO TO 170

  160   CONTINUE
!C  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
!C  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
!C  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN ANY ELEMENT
!C  IN THE SECOND PART, AND THEY MAY NOW BE SORTED INDEPENDENTLY.
  170   COUNT(IFK) = X
        IP = IFK
        GO TO 190

  180   COUNT(I) = X
        IP = I
!C  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  190   IF ((IP-IS).GT. (IF-IP)) GO TO 200
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 210

  200   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
!C  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  210   LNGTH = IF - IS
        IF(LNGTH.LT.0) GO TO 230
        IF(LNGTH.EQ.0) GO TO 260
!C  IF IT CONTAINS MORE THAN ONE ELEMENT STORE IT IN WORKSPACE .
        LA = LA + 2
        MARK(LA) = IF
        MARK(LA-1) = IS
        GO TO 270
!C  IF IT CONTAINS NO ELEMENTS RESELECT THE OTHER SUBDIVISION
!C  AND FIND A DIFFERENT TEST NUMBER. NUMBERS WHICH ARE FOUND TO
!C  EQUAL THE TEST NUMBER ARE SORTED OUT.
  230   IS = MARK(LA-1)
        IF = MARK(LA)
        IFEND = IF - 1
        IK = IS
        DO 240 I = IK,IFEND
          IF (COUNT(I).NE.X) GO TO 250
          IS = IS + 1
  240   CONTINUE
        LA = LA - 2
        GO TO 260

  250   AV = COUNT(I)
        IY = (IF+IS)/2
        COUNT(I) = COUNT(IY)
        COUNT(IY) = AV
        GO TO 270
!C  FIND IF SORTING IS COMPLETED.
  260   IF (LA.LE.0) GO TO 280
!C  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
        IF = MARK(LA)
        IS = MARK(LA-1)
  270 CONTINUE
  280 RETURN

      END







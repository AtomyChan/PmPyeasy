! Subroutine tuneup
! This subroutine takes the parameters in the parameter files and puts them 
! in variables usable by the rest of the subroutines and functions of DoPHOT.
!
! The flags and files are associated with the contents shown below. Notice 
! that the Fortran input and output device units used by DoPHOT follow 
! the same association as the files.
!
!   flag   contents                   file   contents
!    1    PSF type                      1   image_in
!    2    Sky type                      2   image_out
!    3    image_out (y/n)               3   obj_in
!    4    obj_out format                4   obj_out
!    5    obj_in (y/n)                  5   not used - keyboard associated unit
!    6    obj_in format                 6   not used - screen associated unit
!    7    shadow_in (y/n)               7   shadow_in
!    8    shadow_out (y/n)              8   shadow_out
!    9    automatch (y/n)               9   automatch
!   10    automatch format             10   not used
!   11    variable psf (y/n)           11   vpsf coefficients
!   12    ap correc (y/n)              12   apcorr
!   13    ap phot (y/n)                13   logfile
!   14    transformation coeffs (y/n)  14   transformation coeffs
!   15    apply transformations (y/n)  15   finishfile
!   16    obliterate saturated (y/n)   16   not used
!   17    identify volcano stars (y/n) 17   not used

subroutine tuneup(npmax,nff,flags,files,skyguess,fitrad,maskrad,&
     aprad,apskymin,apskymax,eperdn,rnoise,bot,top,tmin,tmax,tfac,&
     nsmax,nit,npstar,nfit0,nfit1,nfit2,nsrmsmin,fac,xpnd,&
     icrit,widobl,cmax,topsat,fsub,fobl,stograt,discrim,sig,chicrit,xtra,&
     crit7,snlim,bumpcrit,sncos,enuff4,enuff7,&
     nbadleft,nbadright,nbadtop,nbadbot,lverb,acc,alim,ava,&
     beta4,beta6,pixthresh,fixpos,fixshape,&
     napcmin,napertures,aperrmax,apmagmaxerr,&
     mapparmag,mapparxy,napsurffit,apfaintmag,&
     npskyp,npskyh,nsskypmin,nsskyhmin,iskymedstep,&
     accskyp,alimskyp,accskyh,alimskyh,itskyp,itskyh,&
     npsffit,npsfmin,ncmin,norder,&
     nsmatchmin,amagmaxm,amagminm,amagerrmaxm,xoffinit,yoffinit,&
     radius,dradius,radmin)

  integer, parameter :: nparam=500
  real,parameter :: fwhmfudge=1.0, rad2deg=57.29578, sig2fw=2.35482
  real :: sig(3),ava(npmax)
  real :: acc(npmax),accskyp(npmax),accskyh(npmax)
  real :: alim(npmax),alimskyp(npmax),alimskyh(npmax)
  character(80) :: prompt,keyword,item,comment,header(nparam)
  character(80) :: flags(nff),files(nff)
  logical autoscale,autoscaleap,autothresh,autoblit,fixpos,fixshape

! Read the parameter values from the default and modified files
  call param_in(nparam,header,nlines)

!   First, read the file names and file usage flags. These variables
!   will be transmitted to the outside world in the array FILES.

!     Read the input image name.  THE NAME MUST EXIST.  
!     If not, it is requested.
  call readitem(nlines,header,'IMAGE_IN',nr,itype,files(1),line)
  ierr=1
  do while (ierr.eq.1)
     if (nr.eq.0) then
        prompt = 'Enter input image name: '
        call query(prompt)
        read(*,'(a)') files(1)
        call readcomment(header(line),comment,ncomm)
        call writeitemc(header(line),'IMAGE_IN',files(1),0,comment,ncomm)
     end if
     call opena(1,files(1),0,ierr)
     close(unit=1)
     if(ierr.eq.1) then
        nr = 0
        print *,'Input image does not exist!'
     end if
  enddo

! Determine if an output picture file is to be saved.
  call readitem(nlines,header,'IMAGE_OUT',nr,itype,files(2),line)
  flags(3)(1:3) = 'NO '
  if(nr.gt.0) then
     flags(3)(1:3) = 'YES'
  end if

! Read input objects list.  Leave null if none is given.
  call readitem(nlines,header,'OBJECTS_IN',nr,itype,files(3),line)
  flags(5)(1:3) = 'NO '
  if(nr.gt.0) then
     flags(5)(1:3) = 'YES'
  end if

! The output objects file name MUST EXIST.  
! If it is not in the parameter file, it is requested.
  call readitem(nlines,header,'OBJECTS_OUT',nr,itype,files(4),line)
  if(nr.eq.0) then
     prompt = 'Enter output objects file name: '
     call query(prompt)
     read(*,'(a)') files(4)
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'OBJECTS_OUT',files(4),0,comment,ncomm)
  end if

! Ask for the shadow file name if it is to be saved.
! If no name is provided, the output file will get a default designation.
  call readitem(nlines,header,'SHADOWFILE_OUT',nr,itype,files(8),line)
  flags(8)(1:3) = 'NO '
  if(nr.gt.0) then
     flags(8)(1:3) = 'YES'
  end if

  flags(7)(1:3) = 'NO '
  if(flags(5)(1:1).eq.'Y') then
     call readitem(nlines,header,'SHADOWFILE_IN',nr,itype,files(7),line)
     if(nr.gt.0) then
        flags(7)(1:3) = 'YES'
     end if
  end if

! Ask for the log file name if it is to be saved.  
! If no name is provided, the output file will get a default designation.
  call readitem(nlines,header,'LOGVERBOSITY',nr,itype,item,line)
  read(item,*) lverb
  if(lverb.gt.0) then
     call readitem(nlines,header,'LOGFILE',nr,itype,files(13),line)
     if(nr.eq.0) then
        files(13) = 'logfile.dat'
        call readcomment(header(line),comment,ncomm)
        call writeitemc(header(line),'LOGFILE',files(13),0,comment,ncomm)
     end if
  end if

! Now read the psf and sky flags.
! If no parameters are present, the defaults shown below are assumed.
  call readitem(nlines,header,'PSFTYPE',nr,itype,flags(1),line)
  if(nr.eq.0) then
     flags(1) = 'PGAUSS'
  else     ! At this time, only PGAUSS has been implemented 
     call upper(flags(1)(1:nr))
     do while (flags(1)(1:5).ne.'PGAUS')
        prompt = 'Enter PSF type: '
        call query(prompt)
        read(*,'(a)') flags(1)
        if(flags(1)(1:5).ne.'PGAUS') then
           print *,'Valid type is PGAUSS.  Try again.'
        end if
     end do
  end if
  call readcomment(header(line),comment,ncomm)
  call writeitemc(header(line),'PSFTYPE',flags(1),0,comment,ncomm)

  call readitem(nlines,header,'SKYTYPE',nr,itype,flags(2),line)
  if(nr.eq.0) then
     flags(2) = 'PLANE'
  else    ! At this time, only PLANE,MEDIAN and HUBBLE have been implemented 
     call upper(flags(2)(1:nr))
     do while ((flags(2)(1:5).ne.'PLANE').and.(flags(2)(1:5).ne.'MEDIA').and.&
          (flags(2)(1:5).ne.'HUBBL'))
        prompt = 'Enter sky type: '
        call query(prompt)
        read(*,'(a)') flags(2)
        if((flags(2)(1:5).ne.'PLANE').and.(flags(2)(1:5).ne.'MEDIA').and.&
          (flags(2)(1:5).ne.'HUBBL')) then
           print *,'Valid types are PLANE, MEDIAN and HUBBLE.  Try again.'
        end if
     end do
  end if
  call readcomment(header(line),comment,ncomm)
  call writeitemc(header(line),'SKYTYPE',flags(2),0,comment,ncomm)

  call readitem(nlines,header,'OBJTYPE_IN',nr,itype,flags(6),line)
  if(nr.eq.0) then
     flags(6) = 'INTERNAL'
  else  ! At this time, only COMPLETE,INTERNAL and BINARY have been implemented
     call upper(flags(6)(1:nr))
     do while ((flags(6)(1:5).ne.'COMPL').and.(flags(6)(1:5).ne.'INTER').and.&
          (flags(6)(1:5).ne.'BINAR'))
        prompt = 'Enter input object file format type: '
        call query(prompt)
        read(*,'(a)') flags(8)
        if((flags(6)(1:5).ne.'COMPL').and.(flags(6)(1:5).ne.'INTER').and.&
             (flags(6)(1:5).ne.'BINAR')) then
           print *,'Valid types are COMPLETE, INTERNAL and BINARY.  Try again.'
        end if
     end do
  end if
  call readcomment(header(line),comment,ncomm)
  call writeitemc(header(line),'OBJTYPE_IN',flags(6),0,comment,ncomm)

  call readitem(nlines,header,'OBJTYPE_OUT',nr,itype,flags(4),line)
  if(nr.eq.0) then
     flags(4) = 'DAOPHOT'
  else
     call upper(flags(4)(1:nr))
     do while ((flags(4)(1:5).ne.'COMPL').and.(flags(4)(1:5).ne.'BINAR').and.&
          (flags(4)(1:5).ne.'INTER').and.(flags(4)(1:5).ne.'INCOM').and.&
          (flags(4)(1:5).ne.'DAOPH'))
        prompt = 'Enter output object file format type: '
        call query(prompt)
        read(*,'(a)') flags(3)
        if  ((flags(4)(1:5).ne.'COMPL').and.(flags(4)(1:5).ne.'BINAR').and.&
             (flags(4)(1:5).ne.'INTER').and.(flags(4)(1:5).ne.'INCOM').and.&
             (flags(4)(1:5).ne.'DAOPH')) then
           print *,'Valid types are COMPLETE, INTERNAL, INCOMPLETE, BINARY, '&
                //'and DAOPHOT. Try again.'
        end if
     end do
  end if
  call readcomment(header(line),comment,ncomm)
  call writeitemc(header(line),'OBJTYPE_OUT',flags(4),0,comment,ncomm)

  call readitem(nlines,header,'FINISHFILE',nr,itype,item,line)
  if(nr.gt.0) read(item,*) files(15)

! The long job now is to extract the necessary keywords one by one. 

  call readitem(nlines,header,'TILT',nr,itype,item,line)
  read(item,*) tilt
  call readitem(nlines,header,'AXIS_RATIO',nr,itype,item,line)
  read(item,*) ar
  call readitem(nlines,header,'FWHM',nr,itype,item,line)
  read(item,*) fwhm
  fwhm = fwhm*fwhmfudge
  tilt = tilt/rad2deg
  sigmasq = (fwhm/sig2fw)**2
  sigfunxsq = sigmasq/(cos(tilt)**2+(sin(tilt)/ar)**2)
  sigfunxy = (-sin(2*tilt)+sin(2*tilt)/ar**2)/(2*sigmasq)
  sigfunysq = sigmasq/((cos(tilt)/ar)**2+sin(tilt)**2)
  fwhmx = sig2fw*sqrt(sigfunxsq)
  fwhmy = sig2fw*sqrt(sigfunysq)
  call readitem(nlines,header,'SKY',nr,itype,item,line)
  read(item,*) skyguess
! Initial parameters for the function used for surface brightness profiles
! See pseud2d.f90 
  ava(1) = skyguess
  ava(2) = 0.0
  ava(3) = 0.0
  ava(4) = 0.0
  ava(5) = sigfunxsq
  ava(6) = sigfunxy
  ava(7) = sigfunysq

! Determine if autoscaling is desired.
  call readitem(nlines,header,'AUTOSCALE',nr,itype,item,line)
  autoscale = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') autoscale = .true.
  if(autoscale) then
     call readitem(nlines,header,'SCALEFITRADIUS',nr,itype,item,line)
     read(item,*) scalefitrad 
     call readitem(nlines,header,'FITRADIUSMIN',nr,itype,item,line)
     read(item,*) fitradmin 
     call readitem(nlines,header,'SCALEMASKRADIUS',nr,itype,item,line)
     read(item,*) scalemaskrad 
     call readitem(nlines,header,'MASKRADIUSMIN',nr,itype,item,line)
     read(item,*) maskradmin
!    Apply these factors to the relevant quantities.
     fitrad = amax1(fwhm*scalefitrad,fitradmin)
     call readitem(nlines,header,'FITRADIUS',nr,itype,item,line)
     write(item,'(f11.2)') fitrad
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'FITRADIUS',item,0,comment,ncomm)
     maskrad = max0(nint(fwhm*scalemaskrad),maskradmin)
     call readitem(nlines,header,'MASKRADIUS',nr,itype,item,line)
     write(item,'(i13)') maskrad
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'MASKRADIUS',item,0,comment,ncomm)
! Else, read the parameters directly from the file.
  else
     call readitem(nlines,header,'FITRADIUS',nr,itype,item,line)
     read(item,*) fitrad
     call readitem(nlines,header,'MASKRADIUS',nr,itype,item,line)
     read(item,*) maskrad
   end if

! Determine if autoscaling of aperture parameters is desired.
  call readitem(nlines,header,'AUTOSCALEAP',nr,itype,item,line)
  autoscaleap = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') autoscaleap = .true.
  if(autoscaleap) then
     call readitem(nlines,header,'SCALEAPRADIUS',nr,itype,item,line)
     read(item,*) scaleaprad 
     call readitem(nlines,header,'APRADIUSMIN',nr,itype,item,line)
     read(item,*) apradmin 
     call readitem(nlines,header,'SCALEAPSKYANNUL',nr,itype,item,line)
     read(item,*) scaleapskyann 
     call readitem(nlines,header,'APSKYANNULUSMIN',nr,itype,item,line)
     read(item,*) apskyannmin
!    Apply these factors to the relevant quantities.
     aprad = amax1(fwhm*scaleaprad,apradmin)
     call readitem(nlines,header,'APRADIUS',nr,itype,item,line)
     write(item,'(f11.2)') aprad
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'APRADIUS',item,0,comment,ncomm)
     call readitem(nlines,header,'APSKYDIST',nr,itype,item,line)
     read(item,*) apskydist
     apskymin = aprad+apskydist
     call readitem(nlines,header,'APSKYMIN',nr,itype,item,line)
     write(item,'(f11.2)') apskymin
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'APSKYMIN',item,0,comment,ncomm)
     apsky = amax1(fwhm*scaleapskyann,apskyannmin)
     apskymax = apskymin+apsky
     call readitem(nlines,header,'APSKYMAX',nr,itype,item,line)
     write(item,'(f11.2)') apskymax
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'APSKYMAX',item,0,comment,ncomm) 
! Else, read the parameters directly from the file.
  else
     call readitem(nlines,header,'APRADIUS',nr,itype,item,line)
     read(item,*) aprad
     call readitem(nlines,header,'APSKYMIN',nr,itype,item,line)
     read(item,*) apskymin
     call readitem(nlines,header,'APSKYMAX',nr,itype,item,line)
     read(item,*) apskymax
   end if

  call readitem(nlines,header,'EPERDN',nr,itype,item,line)
  read(item,*) eperdn
  call readitem(nlines,header,'RDNOISE',nr,itype,item,line)
  read(item,*) rnoise
  call readitem(nlines,header,'TOP',nr,itype,item,line)
  read(item,*) top
  call readitem(nlines,header,'THRESHDEC',nr,itype,item,line)
  read(item,*) tfac
! Ask if auto thresholding is to be done.
  call readitem(nlines,header,'AUTOTHRESH',nr,itype,item,line)
  autothresh = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') autothresh = .true.
  sigskyguess = sqrt(eperdn*skyguess+rnoise**2)/eperdn
  if(autothresh) then
!    Read the needed scaling quantities.
     call readitem(nlines,header,'SIGMABOTTOM',nr,itype,item,line)
     read(item,*) sigbot
     call readitem(nlines,header,'SIGMATHRESHMIN',nr,itype,item,line)
     read(item,*) sigthresh
!    Calculate the sigma of the sky 
     sigskyguess = sqrt(eperdn*skyguess+rnoise**2)/eperdn
     bot = skyguess-sigbot*sigskyguess
     tmin = sigthresh*sigskyguess
     tmax = amax1((top-sky)/2.0,tmin)
!    Apply these parameters.
     call readitem(nlines,header,'BOTTOM',nr,itype,item,line)
     write (item,'(f11.2)') bot
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'BOTTOM',item,0,comment,ncomm)
     call readitem(nlines,header,'THRESHMAX',nr,itype,item,line)
     write (item,'(f11.2)') tmax
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'THRESHMAX',item,0,comment,ncomm)
     call readitem(nlines,header,'THRESHMIN',nr,itype,item,line)
     write (item,'(f11.2)') tmin
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'THRESHMIN',item,0,comment,ncomm)
! Else, read the parameters directly from the file.
  else
     call readitem(nlines,header,'BOTTOM',nr,itype,item,line)
     read(item,*) bot
     call readitem(nlines,header,'THRESHMAX',nr,itype,item,line)
     read(item,*) tmax  
     call readitem(nlines,header,'THRESHMIN',nr,itype,item,line)
     read(item,*) tmin
  end if

  call readitem(nlines,header,'FSUB',nr,itype,item,line)
  read(item,*) fsub
  call readitem(nlines,header,'FOBL',nr,itype,item,line)
  read(item,*) fobl
! Ask if auto obliteration is to be done.
  call readitem(nlines,header,'AUTOBLIT',nr,itype,item,line)
  autoblit = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') autoblit = .true.
  if(autoblit) then
     call readitem(nlines,header,'ICRITMIN',nr,itype,item,line)
     read(item,*) icritmin
     call readitem(nlines,header,'SCALEICRIT',nr,itype,item,line)
     read(item,*) scaleicrit
     call readitem(nlines,header,'SCALECENT',nr,itype,item,line)
     read(item,*) scalecent
     call readitem(nlines,header,'SCALETOPSAT',nr,itype,item,line)
     read(item,*) scaletopsat
     icrit = max0(icritmin,nint(scaleicrit*((0.5*fwhm)**2)))
     cmax = scalecent*top
     topsat = scaletopsat*top

     call readitem(nlines,header,'ICRIT',nr,itype,item,line)
     write(item,'(i8)') icrit
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'ICRIT',item,0,comment,ncomm)
     call readitem(nlines,header,'CENTINTMAX',nr,itype,item,line)
     write(item,'(f11.2)') cmax
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'CENTINTMAX',item,0,comment,ncomm)
     call readitem(nlines,header,'TOPSAT',nr,itype,item,line)
     write(item,'(f11.2)') topsat
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'TOPSAT',item,0,comment,ncomm)
  else
     call readitem(nlines,header,'ICRIT',nr,itype,item,line)
     read(item,*) icrit
     call readitem(nlines,header,'CENTINTMAX',nr,itype,item,line)
     read(item,*) cmax
     call readitem(nlines,header,'TOPSAT',nr,itype,item,line)
     read(item,*) topsat
  end if

! Read the variables used to describe the sky functions
  call readitem(nlines,header,'NPSKYP',nr,itype,item,line)
  read (item,*) npskyp
  call readitem(nlines,header,'NPSKYH',nr,itype,item,line)
  read (item,*) npskyh
  call readitem(nlines,header,'NSSKYPMIN',nr,itype,item,line)
  read (item,*) nsskypmin
  call readitem(nlines,header,'NSSKYHMIN',nr,itype,item,line)
  read (item,*) nsskyhmin
  call readitem(nlines,header,'ITSKYP',nr,itype,item,line)
  read (item,*) itskyp
  call readitem(nlines,header,'ITSKYH',nr,itype,item,line)
  read (item,*) itskyh
  call readitem(nlines,header,'MEDSTEP',nr,itype,item,line)
  read(item,*) iskymedstep 
  call readitem(nlines,header,'AUTOMEDSCALE',nr,itype,item,line)
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') then
     call readitem(nlines,header,'SCALEMEDSTEP',nr,itype,item,line)
     read(item,*) scalemedstep 
     iskymedstep = max0(int(fwhm*scalemedstep),iskymedstep)
  end if

!  Ask if positions and shapes will be fixed.
  call readitem(nlines,header,'FIXPOS',nr,itype,item,line)
  fixpos = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') fixpos = .true.
  call readitem(nlines,header,'FIXSHAPE',nr,itype,item,line)
  fixshape = .false.
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') fixshape = .true.

!  Variable PSF?
  call readitem(nlines,header,'VARIABLE_PSF',nr,itype,item,line)
  flags(11) = 'NO'
  npsffit = 1
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') flags(11) = 'YES'
  if(flags(11).eq.'YES') then
     call readitem(nlines,header,'NPSFORDER',nr,itype,item,line)
     read(item,*) ipsford
     if(ipsford.eq.1) npsffit = 3
     if(ipsford.eq.2) npsffit = 6
     if(ipsford.eq.3) npsffit = 10
     if(ipsford.eq.4) npsffit = 15
     if(ipsford.eq.5) npsffit = 21
     if(ipsford.eq.6) npsffit = 28
     if(ipsford.eq.7) npsffit = 36
     if(ipsford.eq.8) npsffit = 45
     if(ipsford.eq.9) npsffit = 55
     if(ipsford.ge.10) npsffit = 66
     call readitem(nlines,header,'NPSFMIN',nr,itype,item,line)
     read(item,*) npsfmin 
  end if
  if(ipsford.eq.0) flags(11) = 'NO'
  call readitem(nlines,header,'LOGVPSF',nr,itype,item,line)
  if(nr.gt.0) read (item,'(a)') files(11)

! Just aperture photometry?
  call readitem(nlines,header,'APERTURE_PHOT',nr,itype,item,line)
  flags(13) = 'NO'
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') flags(13) = 'YES'

! Mask/obliterate saturated stars?
  call readitem(nlines,header,'MASK_SATURATED',nr,itype,item,line)
  flags(16) = 'YES'
  if(item(1:1).eq.'n'.or.item(1:1).eq.'N') flags(16) = 'NO'

! Better identify 'volcano' stars?
  call readitem(nlines,header,'VOLCANO',nr,itype,item,line)
  flags(17) = 'NO'
  if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') flags(17) = 'YES'

! Matching files?
  flags(9) = 'NO'
  call readitem(nlines,header,'AUTOMATCH_FILE',nr,itype,item,line)
  if(nr.gt.0) then
     read(item,*) files(9)
     flags(9) = 'YES'
  end if
  flags(15) = 'NO'
  if(flags(9)(1:1).eq.'Y') then
     call readitem(nlines,header,'APPLY_TRANS',nr,itype,item,line)
     if(item(1:1).eq.'y'.or.item(1:1).eq.'Y') flags(15) = 'YES'
     call readitem(nlines,header,'AUTOMATCH_TYPE',nr,itype,item,line)
     if(nr.eq.0) then
        flags(10) = 'INTERNAL'
     else
        call upper(item(1:nr))
        read(item,*) flags(10)
        do while ((flags(10)(1:5).ne.'COMPL').and.(flags(10)(1:5).ne.'INTER')&
             .and.(flags(10)(1:5).ne.'BINAR').and.(flags(10)(1:5).ne.'DAOPH'))
           prompt = 'Enter automatch file format type: '
           call query(prompt)
           read(*,'(a)') flags(10)
           if ((flags(10)(1:5).ne.'COMPL').and.(flags(10)(1:5).ne.'INTER')&
                .and.(flags(10)(1:5).ne.'BINAR')&
                .and.(flags(10)(1:5).ne.'DAOPH')) then
              print *,'Valid types are COMPLETE, INTERNAL, BINARY and DAOPHOT.'&
                   //'Try again.'
           end if
        end do
     end if
     call readcomment(header(line),comment,ncomm)
     call writeitemc(header(line),'AUTOMATCH_TYPE',flags(10),0,&
          comment,ncomm)
  end if



  call readitem(nlines,header,'RESIDNOISE',nr,itype,item,line)
  read(item,*) fac
  
  call readitem(nlines,header,'FOOTPRINT_NOISE',nr,itype,item,line)
  read(item,*) xpnd

  call readitem(nlines,header,'STARGALKNOB',nr,itype,item,line)
  read(item,*) stograt
  call readitem(nlines,header,'STARCOSKNOB',nr,itype,item,line)
  read(item,*) discrim
  
  call readitem(nlines,header,'SNLIM7',nr,itype,item,line)
  read(item,*) crit7
  call readitem(nlines,header,'SNLIM',nr,itype,item,line)
  read(item,*) snlim
  call readitem(nlines,header,'SNLIMMASK',nr,itype,item,line)
  read(item,*) bumpcrit
  call readitem(nlines,header,'SNLIMCOS',nr,itype,item,line)
  read(item,*) sncos
  
  call readitem(nlines,header,'NBADLEFT',nr,itype,item,line)
  read(item,*) nbadleft
  call readitem(nlines,header,'NBADRIGHT',nr,itype,item,line)
  read(item,*) nbadright
  call readitem(nlines,header,'NBADTOP',nr,itype,item,line)
  read(item,*) nbadtop
  call readitem(nlines,header,'NBADBOT',nr,itype,item,line)
  read(item,*) nbadbot

  call readitem(nlines,header,'NCMIN',nr,itype,item,line)
  read(item,*) ncmin 
  if(ncmin.eq.0) ncmin = 2

  call readitem(nlines,header,'NSMAX',nr,itype,item,line)
  read(item,*) nsmax
  call readitem(nlines,header,'NFITITER',nr,itype,item,line)
  read(item,*) nit
  call readitem(nlines,header,'NPARAM',nr,itype,item,line)
  read(item,*) npstar
  if (npstar.gt.npmax) then
     print*,'Change npstar in dophot or in .tuneup/paramdefault'
     stop
  endif
  call readitem(nlines,header,'NFITMAG',nr,itype,item,line)
  read(item,*) nfit0
  call readitem(nlines,header,'NFITPOS',nr,itype,item,line)
  read(item,*) nfit1
  call readitem(nlines,header,'NFITSHAPE',nr,itype,item,line)
  read(item,*) nfit2
  call readitem(nlines,header,'NSRMSMIN',nr,itype,item,line)
  read (item,*) nsrmsmin
  call readitem(nlines,header,'CHI2MINBIG',nr,itype,item,line)
  read(item,*) chicrit
  call readitem(nlines,header,'XTRA',nr,itype,item,line)
  read(item,*) xtra
  call readitem(nlines,header,'SIGMA1',nr,itype,item,line)
  read(item,*) sig(1)
  call readitem(nlines,header,'SIGMA2',nr,itype,item,line)
  read(item,*) sig(2)
  call readitem(nlines,header,'SIGMA3',nr,itype,item,line)
  read(item,*) sig(3)
  call readitem(nlines,header,'ENUFF4',nr,itype,item,line)
  read(item,*) enuff4
  call readitem(nlines,header,'ENUFF7',nr,itype,item,line)
  read(item,*) enuff7
  call readitem(nlines,header,'COSOBLSIZE',nr,itype,item,line)
  read(item,*) widobl
  call readitem(nlines,header,'PIXTHRESH',nr,itype,item,line)
  read(item,*) pixthresh
  call readitem(nlines,header,'BETA4',nr,itype,item,line)
  read(item,*) beta4
  call readitem(nlines,header,'BETA6',nr,itype,item,line)
  read(item,*) beta6

! Ask for aperture correction parameters  
  flags(12)(1:3) = 'NO '
  call readitem(nlines,header,'APCORRFILE',nr,itype,item,line)
  read (item,'(a)') files(12)
  flags(12)(1:3) = 'YES'
  call readitem(nlines,header,'NAPCMIN',nr,itype,item,line)
  read (item,*) napcmin
  call readitem(nlines,header,'NAPERTURES',nr,itype,item,line)
  read (item,*) napertures
  call readitem(nlines,header,'APERRMAX',nr,itype,item,line)
  read(item,*) aperrmax
  call readitem(nlines,header,'MAPPARMAG',nr,itype,item,line)
  read (item,*) mapparmag
  call readitem(nlines,header,'MAPPARXY',nr,itype,item,line)
  read (item,*) mapparxy
  call readitem(nlines,header,'NAPSURFFIT',nr,itype,item,line)
  read (item,*) napsurffit
  call readitem(nlines,header,'APFAINTMAG',nr,itype,item,line)
  read (item,*) apfaintmag
  call readitem(nlines,header,'APMAG_MAXERR',nr,itype,item,line)
  read(item,*) apmagmaxerr
     
! Ask for automatching parameters
  call readitem(nlines,header,'NSMATCHMIN',nr,itype,item,line)
  read(item,*) nsmatchmin
  call readitem(nlines,header,'AMAGMIN_MATCH',nr,itype,item,line)
  read(item,*) amagminm
  call readitem(nlines,header,'AMAGMAX_MATCH',nr,itype,item,line)
  read(item,*) amagmaxm
  if(amagminm.gt.amagmaxm) then
     saco = amagminm
     amagminm = amagmaxm
     amagmaxm = saco
  end if
  call readitem(nlines,header,'AMAGERRMAX_MATCH',nr,itype,item,line)
  read(item,*) amagerrmaxm
  call readitem(nlines,header,'STARTRAD',nr,itype,item,line)
  read(item,*) radius
  call readitem(nlines,header,'DELTARAD',nr,itype,item,line)
  read(item,*) dradius
  call readitem(nlines,header,'FINALRAD',nr,itype,item,line)
  read(item,*) radmin
  call readitem(nlines,header,'XOFFSET',nr,itype,item,line)
  read(item,*) xoffinit
  call readitem(nlines,header,'YOFFSET',nr,itype,item,line)
  read(item,*) yoffinit
  call readitem(nlines,header,'NTRANSORDER',nr,itype,item,line)
  read(item,*) norder
  call readitem(nlines,header,'TRANS_COEFS',nr,itype,files(14),line)
  if(nr.gt.0) then
     flags(14) = 'YES'
  else
     flags(14) = 'NO'
  end if

! Parameter limits   
  keyword(1:6) = 'RELACC'
  do i=1,npstar
     write(keyword(7:7),'(i1)') i
     call readitem(nlines,header,keyword(1:7),nr,itype,item,line)
     read(item,*) acc(i)
  end do
 
  keyword(1:6) = 'ABSLIM'
  do i=1,npstar
     write(keyword(7:7),'(i1)') i
     call readitem(nlines,header,keyword(1:7),nr,itype,item,line)
     read(item,*) alim(i)
  end do
  
  keyword(1:14) = 'RELACCSKYPLANE'
  do i=1,npskyp
     write(keyword(15:15),'(i1)') i
     call readitem(nlines,header,keyword(1:15),nr,itype,item,line)
     read(item,*) accskyp(i)
  end do
 
  keyword(1:14) = 'ABSLIMSKYPLANE'
  do i=1,npskyp
     write(keyword(15:15),'(i1)') i
     call readitem(nlines,header,keyword(1:15),nr,itype,item,line)
     read(item,*) alimskyp(i)
  end do
  
  keyword(1:12) = 'RELACCSKYHUB'
  do i=1,npskyh
     write(keyword(13:13),'(i1)') i
     call readitem(nlines,header,keyword(1:13),nr,itype,item,line)
     read(item,*) accskyh(i)
  end do
 
  keyword(1:12) = 'ABSLIMSKYH'
  do i=1,npskyh
     write(keyword(13:13),'(i1)') i
     call readitem(nlines,header,keyword(1:13),nr,itype,item,line)
     read(item,*) alimskyh(i)
  end do

  

! Write current parameters in a file.
  call param_out(header,nlines)
    
end subroutine tuneup



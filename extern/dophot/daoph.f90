! This file contains the subroutines to write files 
! with the DAOPHOT dophot format
!
! DAOPH_OUT
! DAOPH_HEAD
! DAOPH_HEADINIT
! DAOPH_FINAL
! DAOPH_IN
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine daoph_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)
! This subroutine writes the data with DAOPHOT format.

subroutine daoph_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)

  real starpar(npstar),appar(npap)
  logical applyap
  character(120) :: fmt

  if ((itype.eq.4).or.(itype.eq.6).or.(itype.eq.8)) then
     fmag = 99.999
     error = 9.999
     chi = 0.0
     apcorr = 0.0
  else
     if(starpar(2).le.0.0) return
     fmag=SHAPEMAG(starpar,npstar)
     if(fmag.lt.50.0) then
        error = appar(4)
        apcorr = appar(3)
        if(applyap) fmag = fmag+apcorr
!       A quasi chi-sq like statistic of some usefulness.
!       chi=error/sigmae 
!       where error=appar(4) and sigmae=1/sqrt(3.0*starpar(2))
        chi = appar(4)*sqrt(3.0*starpar(2))
     else
        error = 9.999
        chi = 0.0
        apcorr = 0.0
     end if
  end if

! Order is: No.,xpos,ypos,fitmag,err_fitmag,
!           fitsky,obstype,chi,apcorr(diff_fit_ap)
  fmt="(i6,2f9.2,2f9.3,f9.2,f9.0,f9.2,f9.3)"
  write(iunit,fmt) i,starpar(3),starpar(4),fmag,error,&
       starpar(1),float(itype),chi,apcorr

end subroutine daoph_out

!---------------------------------------------------------------------------
!
! Subroutine daoph_head(iunit,census,ncensus,ntot,fwhm,sky)
! This subroutine writes the header of a file with DAOPHOT format.

subroutine daoph_head(iunit,census,ncensus,ntot,fwhm,sky)

  integer :: census(ncensus)
  character(2) :: fmtvar
  character(120) :: fmt1,fmt2
  
  rewind(iunit)

  write(fmtvar,'(i2)') ncensus
  fmt1="(a4,"//fmtvar//"i6,a7,a6,a9)"
  fmt2="(4x,"//fmtvar//"i6,i7,f6.1,f9.1)"
  write(iunit,fmt1) ' NL ',(k,k=1,ncensus),'Total','FWHM','Sky'
  write(iunit,fmt2) (census(k),k=1,ncensus),ntot,fwhm,sky
  write(iunit,*) ' '

end subroutine daoph_head

!---------------------------------------------------------------------------
!
! Subroutine daoph_headinit(iunit,ntype)
! This subroutine writes an initial header for a file with DAOPHOT format.

subroutine daoph_headinit(iunit,ntype)

  integer :: census(ntype)

  ntot = 0
  census = 0
  fwhm = 0.0
  sky = 0.0
  call daoph_head(iunit,census,ntype,ntot,fwhm,sky)

end subroutine daoph_headinit

!---------------------------------------------------------------------------
!
! Subroutine daoph_final(iunit,ntype,fileout,nstot,npstar,ava)
! This subrotine writes the final (and correct) header for 
! a file with DAOPHOT format, and orders the data according to the y coordinate
 
subroutine daoph_final(iunit,ntype,fileout,nstot,npstar,ava)

  character(80) :: lines(nstot),fileout
  real :: y(nstot),ava(npstar)
  integer :: iindex(nstot),census(ntype)

  census=0

  call opena(iunit,fileout,2,ierr)
! Get rid of the header
  do i=1,3
     read(iunit,*)
  end do
! Read and count the data.
  do i=1,nstot
     read(iunit,'(a)',err=215,end=215) lines(i)
     read(lines(i)(16:24),'(f9.2)') y(i)
     read(lines(i)(52:60),'(f9.0)') typpe
     itype = int(typpe)
     j = mod(itype,10)
     census(j) = census(j)+1
  end do
215 ntot = i-1

! Write the final header correctly
  call ellipse(ava(5),ava(6),ava(7),dum,fwhmx,fwhmy,dum)
  fwhm = (fwhmx+fwhmy)/2.0
  sky = ava(1)
  call daoph_head(iunit,census,ntype,ntot,fwhm,sky)
! Sort the data by y and write them on the file
  call quick(y,ntot,iindex)
  do i=1,ntot
     ii = iindex(i)
     read(lines(ii)(52:60),'(f9.0)') typpe
     itype = int(typpe)
     j = mod(itype,10)
     if((j.eq.4).or.(j.eq.6).or.(j.eq.8)) cycle
     read(lines(ii)(25:42),'(2f9.3)') fmag,error
     if((fmag.gt.50.0).or.(error.le.0.0)) cycle
!     if((j.eq.9).and.(error.gt.0.05)) cycle
     write(lines(ii)(1:6),'(i6)') i
     write(iunit,'(a)') lines(ii)
  end do
  close (unit=iunit)

end subroutine daoph_final

!---------------------------------------------------------------------------
!
! Subroutine daoph_in(i,itype,starpar,npstar,appar,npap,iunit,applyap)
! This subroutine reads a file with DAOPHOT format.
! Notice that the information recovered is very restricted
! and therefore it should only be use for automatching.
! Notice also that we should get rid of the header before using this subroutine
 
subroutine daoph_in(i,itype,starpar,npstar,appar,npap,fmag,ios,iunit)

  real :: starpar(npstar),appar(npap)

  read(iunit,*,iostat=ios) i,starpar(3),starpar(4),fmag,appar(4),&
       starpar(1),rtype,chi,apcorr

  itype=nint(rtype)

end subroutine daoph_in


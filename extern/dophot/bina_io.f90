! This file contains the subroutines to read or write the BINARY dophot format
!
! BINA_IN
! BINA_OUT
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine bina_in(i,itype,starpar,npstar,appar,npap,ios,iunit)
! This subroutine reads a file with BINARY format.

subroutine bina_in(i,itype,starpar,npstar,appar,npap,ios,iunit)

  real :: starpar(npstar),appar(npap)

  write(iunit,rec=i,iostat=ios) itype,xc,yc,fmag,appar(4),starpar(1),&
       amajor,aminor,tilt
  if(ios.ne.0) return

! Convert position centers so that the center of the 1st pixel is 1.0 
! (as required in the internal representation) rather than
! 0.5 as used in the COMPLETE style output
  starpar(3) = xc+0.5
  starpar(4) = yc+0.5

  if(itype.eq.8) then
     if((tilt.gt.89.5).and.(tilt.lt.90.5)) then
        starpar(5) = aminor
        starpar(7) = amajor
        starpar(6) = 0.0
        starpar(2) = 0.0
     else
        starpar(5) = amajor
        starpar(7) = aminor
        starpar(6) = 0.0
        starpar(2) = 0.0
     end if
  else
     call ellipint(amajor,aminor,tilt,area,starpar(5),starpar(6),starpar(7))
     if(fmag.gt.50.0) then
        starpar(2) = 0.0
     else
        starpar(2) = 10.**(-0.4*(fmag-23.5))
        starpar(2) = starpar(2)/area
     end if
  end if

end subroutine bina_in

!---------------------------------------------------------------------------
!
! Subroutine bina_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)
! This subroutine writes a file with BINARY format.

subroutine bina_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)

  real :: starpar(npstar), appar(npap)
  logical :: applyap

  if(itype.ne.8) then
     call ellipse(starpar(5),starpar(6),starpar(7),area,amajor,aminor,tilt)
     fmag = area*starpar(2)
     if(fmag.le.0.0) then
        fmag = 99.999
     else 
        fmag = -2.5*log10(fmag)
     end if
  else
     if(starpar(6).eq.-1.0) then
        fmag = 99.999
     else  
        fmag = -99.999
     end if
     if(starpar(5).ge.starpar(7)) then
        amajor = starpar(5)
        aminor = starpar(7)
        tilt = 0.0
     else 
        amajor = starpar(7)
        aminor = starpar(5)
        tilt = 90.0
     end if
  end if
  fmag = fmag + 23.5
  
  if(applyap.and.(abs(fmag).lt.90.0)) fmag = fmag + appar(3)

! Fix the co-ordinates so that the center of the first pixel is 0.5
! and not 1.0 as in the internal representation
  xc = starpar(3)-0.5
  yc = starpar(4)-0.5

! Order is: obtype,xpos,ypos,fitmag,err_fitmag,fitsky,FWHM_major,FWHM_minor,tilt
  write(iunit,rec=i) itype,xc,yc,fmag,appar(4),starpar(1),amajor,aminor,tilt

end subroutine bina_out

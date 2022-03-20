! Subroutine incomp_out(i,itype,starpar,npstar,appar,npap,probg,iunit,applyap)
! This subroutine writes a file with INCOMPLETE format.

subroutine incomp_out(i,itype,starpar,npstar,appar,npap,probg,iunit,applyap)

  real :: starpar(npstar),appar(npap)
  logical :: applyap
  character(120) :: fmt

  if(itype.ne.8) then
     area = ELLIPAREA(starpar(5),starpar(6),starpar(7))
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
  end if

  if(applyap.and.(abs(fmag).lt.90.0)) fmag = fmag + appar(3)

! Fix the co-ordinates so that the center of the first pixel is 0.5
! and not 1.0 as in the internal representation
  xc = starpar(3)-0.5
  yc = starpar(4)-0.5

!Order is:id,xpos,ypos,fitmag,err_fitmagfitsky,obstype,probg,apcorr(diff_fit_ap)
  fmt="(i6,2f9.2,2f9.3,f9.2,f9.0,f9.2,f9.3)"
  write(iunit,fmt) i,xc,yc,fmag,appar(4),starpar(1),float(itype),probg,appar(3)
   
end subroutine incomp_out


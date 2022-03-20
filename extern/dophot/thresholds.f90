! Subroutine thresholds(tmin,tmax,tfac,mthresh > thresh,ntsh)
! This subroutine calculates the different thresholds that we are gonna use
! to find stars.

subroutine thresholds(tmin,tmax,tfac,mthresh,thresh,ntsh)

  real :: thresh(mthresh),hold(mthresh)

  i = 1
  do
     hold(i) = tmin*2.0**(tfac*float(i-1))
     if(hold(i).gt.tmax) exit
     i = i + 1
  end do
  ntsh = i - 1

  if(ntsh.eq.0) then
     print *,'TMAX < TMIN; Exiting.'
     call exit(2)       ! subroutine exit is an extension to f90
  end if

  do i=1,ntsh
     ii = ntsh + 1 - i
     thresh(ii) = hold(i)
  end do
  
end subroutine thresholds

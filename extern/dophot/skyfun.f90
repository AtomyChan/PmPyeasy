! This file contains the subroutine and functions that describe 
! the differenet sky models:
!
! SKYFUN_PLANE        Subroutine
! SKYFUN_HUB          Subroutine
! SKYFUN_MED          Function
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine skyfun_plane(xy,a,npar,dum1,dum2, > fa,skyplaneval)
!
! This subroutine models the variation accros the sky using 
! a uniform gradient model.

subroutine skyfun_plane(xy,a,npar,dum1,dum2,fa,skyplaneval)
 
  real :: xy(2)
  real :: a(npar),fa(npar)
  real :: dum1,dum2     ! Only neccessary for coherence in subrout. chisq
    
  fa(1) = 1
  fa(2) = xy(1)
  fa(3) = xy(2)

  skyplaneval = 0
  do i = 1, 3
     skyplaneval = skyplaneval+a(i)*fa(i)
  end do

end subroutine skyfun_plane

!---------------------------------------------------------------------------
!
! Subroutine skyfun_hub(xy,a,npar,dum1,dum2, > fa,skyhubval)
!
! This subroutine models the variation accros the sky using 
! a modified Hubble profile.

subroutine skyfun_hub(xy,a,npar,dum1,dum2,fa,skyhubval)

  real, parameter :: half=0.5,expmin=-23. 
  real :: xy(2)
  real :: a(npar),fa(npar)
  real :: dum1,dum2     ! Only neccessary for coherence in subrout. chisq

  x = xy(1)-a(3)
  y = xy(2)-a(4)
  a5 = 1./a(5)
  a7 = 1./a(7)
  t5 = a5*x
  t6 = a(6)*y
  t7 = a7*y
  t1 = half*((t5+2*t6)*x+t7*y)
  if (t1.gt.0) then
     den = 1.+t1
     ddendt1 = 1.               !Partial derivative
     pexp = 1./den
  else
     t1 = amax1(t1,expmin)
     pexp = exp(-t1)
     den = 1.
     ddendt1 = 1.
  end if

  skyhubval = a(2)*pexp+a(1)
  
  if(den.eq.0) write(6,*) 'Divide by Zero 18'
  dIdt1 = a(2)*pexp*ddendt1/den     !Actually -dIdt1, but taken care in fa()
  ! fa() is the derivative of function I(x,y) with respect to the parameters
  fa(1) = 1.0
  fa(2) = pexp
  fa(3) = (t5+t6)*dIdt1
  fa(4) = (a(6)*x+t7)*dIdt1
  fa(5) = half*t5**2*dIdt1
  fa(6) = -x*y*dIdt1
  fa(7) = half*t7**2*dIdt1

end subroutine skyfun_hub

!-----------------------------------------------------------------------
!
! Function skyfun_med(xy,a,nfast,nslow,iskymedstep)
! This subroutine models the variation accros the sky using 
! a median filter.

function skyfun_med(xy,a,nfast,nslow,iskymedstep)

  real ::  xy(2),a(nfast/iskymedstep,nslow/iskymedstep)

  mx = nfast/iskymedstep
  my = nslow/iskymedstep
  xskymed = xy(1)/iskymedstep
  yskymed = xy(2)/iskymedstep
  ix = int(xskymed)
  iy = int(yskymed)


  if(xskymed.ge.1.0.and.xskymed.lt.float(mx).and.&
       yskymed.ge.1.0.and.yskymed.lt.float(my)) then
     r1 = (xskymed-ix)*(a(ix+1,iy)-a(ix,iy))+a(ix,iy)
     r2 = (xskymed-ix)*(a(ix+1,iy+1)-a(ix,iy+1))+a(ix,iy+1)
     skyfun_med = ((yskymed-iy)*(r2-r1)) + r1
  else if(xskymed.lt.1.0) then
     if(yskymed.lt.1.0) then
        skyfun_med = a(1,1)
     else if(yskymed.ge.float(my)) then
        skyfun_med = a(1,my)
     else
        skyfun_med=(yskymed-iy)*(a(1,iy+1)-a(1,iy))+a(1,iy)
     end if
  else if(xskymed.ge.float(mx)) then
     if(yskymed.lt.1.0) then
        skyfun_med = a(mx,1)
     else if(yskymed.ge.float(my)) then
        skyfun_med = a(mx,my)
     else
        skyfun_med=(yskymed-iy)*(a(mx,iy+1)-a(mx,iy))+a(mx,iy)
     end if
  else if(yskymed.lt.1.0.and.xskymed.ge.1.0.and.xskymed.lt.float(mx)) then
     skyfun_med=(xskymed-ix)*(a(ix+1,1)-a(ix,1))+a(ix,1)
  else if(yskymed.ge.float(my).and.xskymed.ge.1.0.and.xskymed.lt.float(mx)) then
     skyfun_med=(xskymed-ix)*(a(ix+1,my)-a(ix,my))+a(ix,my)
  end if

end function skyfun_med

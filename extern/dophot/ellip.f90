! This file contains:
!
! ELLIPSE       Subroutine
! ELLIPINT      Subroutine
! ELLIPAREA     Function
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine ellipse(a5,a6,a7, > area,amajor,aminor,angle)
! This subroutine provides the characteristics (semiaxis, area, and tilt)
! of the elliptical footprint of an object from its 'shape' parameters

subroutine ellipse(a5,a6,a7,area,amajor,aminor,angle)

  real,parameter :: pi=3.14159, rad2deg=57.29578, sig2fw=2.35482

! 1/2*(x2/A+2Bxy+y2/C)=1 ==> Ax2+Bxy+Cy2=1 
  a = 1.0/(2.0*a5)
  b = a6
  c = 1.0/(2.0*a7)

  angle = rad2deg*atan2(b,c-a)/2.0    !See ellipint

! Ax2+Bxy+Cy2=1 ==> Ru2+Sv2=1
  root = sqrt((a-c)**2+b**2)
  root1 = (a+c+root)/2.0
  root2 = root1-root          !root2=(a+c-root)/2

! Remember that the area of a ellipse of the form x2/p2+y2/q2=1 is pi*p*q
  if (root2.gt.0.0) then
     area = pi/sqrt(root1*root2)
  else
     area = 0.0
  end if

  amajor = sig2fw*sqrt(1./(2.0*root2))
  aminor = sig2fw*sqrt(1./(2.0*root1))

end subroutine ellipse

!---------------------------------------------------------------------------
!
! Subroutine ellipint(amajor,aminor,tilt, > area,a5,a6,a7)
! Inverse of the subroutine ellipse (see above)

subroutine ellipint(amajor,aminor,tilt,area,a5,a6,a7)

  real,parameter :: pi=3.14159, rad2deg=57.29578, sig2fw=2.35482

  if(amajor.lt.aminor) then 
     print *,'amajor,aminor = ',amajor,aminor
     print *,' maj-axis .lt. min-axis inverted in input data'
     test=amajor
     amajor=aminor
     aminor=test
  end if

  root1 = ((sig2fw/aminor)**2)/2
  root2 = ((sig2fw/amajor)**2)/2

  area = pi/sqrt(root1*root2)

  b = (root1-root2)*sin(2.*(tilt/rad2deg))
  c_mn_a = (root1-root2)*cos(2.*(tilt/rad2deg))      !c-a
  c_pl_a = root1+root2                               !c+a
  c = 0.5*(c_pl_a+c_mn_a)
  a = 0.5*(c_pl_a-c_mn_a)

! 1/2(x2/A+2Bxy+y2/C)=1 <== Ax2+Bxy+Cy2=1 
  a5 = 1.0/(2.0*a)
  a6 = b
  a7 = 1.0/(2.0*c)

end subroutine ellipint

!---------------------------------------------------------------------------
!
! Function elliparea(a5,a6,a7)
! This function provides the area of the elliptical footprint of an object 
! from its 'shape' parameters. The eq. used in subroutine ellipse, 
! after some transformations, becomes the one shown here. 
! Since this is faster, it is the one we use.

function elliparea(a5,a6,a7)

  real,parameter :: pi=3.14159

! 1/2(x2/A+2Bxy+y2/C)=1 ==> Ax2+Bxy+Cy2=1 
  a = 1.0/(2.0*a5)
  b = a6
  c = 1.0/(2.0*a7)

  delta = 4*a*c-b**2
  if (delta.gt.0) then
     elliparea = 2*pi/sqrt(delta)
  else
     elliparea = 0.0
  end if
  
end function elliparea

! This file contains:
!
! SURFACE           Subroutine
! SURFACE1          Subroutine
! SUFPOLY           Subroutine
! SURFEVAL          Function
! TERM              Function
! GAUSSJ            Subroutine
! MGETMFIT          Function
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine surface(xy,aparam,mfit,> fparam,surfaceval)
! This subroutine calculates the values for a polynomial bivariate function.

subroutine surface(xy,aparam,mfit,par1,par2,fparam,surfaceval)

  integer, parameter ::  mmax=66
  real :: xy(2)
  real :: aparam(mfit),fparam(mfit)
  real :: par1,par2     !Only neccessary for coherence in subrout. chisq
 
  if (mfit.gt.mmax) then
     print *,'mfit >  mmax in subrout.surface. Exiting'
    call exit(5)    ! Subrout exit is an extension to f90
  end if
    
  x = xy(1)
  y = xy(2)

  surfaceval = aparam(1)
  if(mfit.ge.3) surfaceval = surfaceval + aparam(2)*x + aparam(3)*y 
  if(mfit.ge.6) then
     x2 = x*x
     y2 = y*y
     surfaceval = surfaceval + aparam(4)*x2 + aparam(5)*x*y + aparam(6)*y2
  end if
  if(mfit.ge.10) then
     x3 = x2*x
     x2y = x2*y
     xy2 = x*y2
     y3 = y2*y
     surfaceval = surfaceval + aparam(7)*x3 + aparam(8)*x2y + &
          aparam(9)*xy2 + aparam(10)*y3
  end if
  if(mfit.ge.15) then
     x4 = x3*x
     x3y = x3*y
     x2y2 = x2*y2
     xy3 = x*y3
     y4 = y*y3
     surfaceval = surfaceval + aparam(11)*x4 + aparam(12)*x3y + &
          aparam(13)*x2y2 + aparam(14)*xy3 + aparam(15)*y4
  end if
  if(mfit.ge.21) then
     x5 = x4*x
     x4y = x4*y
     x3y2 = x3*y2
     x2y3 = x2*y3
     xy4 = x*y4
     y5 = y*y4
     surfaceval = surfaceval + aparam(16)*x5 + aparam(17)*x4y + &
          aparam(18)*x3y2 + aparam(19)*x2y3 + aparam(20)*xy4 + aparam(21)*y5
  end if
  if(mfit.ge.28) then
     x6 = x5*x
     x5y = x5*y
     x4y2 = x4*y2
     x3y3 = x3*y3
     x2y4 = x2*y4
     xy5 = x*y5
     y6 = y*y5
     surfaceval = surfaceval+aparam(22)*x6+ aparam(23)*x5y + aparam(24)*x4y2 +&
          aparam(25)*x3y3 + aparam(26)*x2y4 + aparam(27)*xy5 + aparam(28)*y6
  end if
  if(mfit.ge.36) then
     x7 = x6*x
     x6y = x6*y
     x5y2 = x5*y2
     x4y3 = x4*y3
     x3y4 = x3*y4
     x2y5 = x2*y5
     xy6 = x*y6
     y7 = y6*y
     surfaceval = surfaceval + aparam(29)*x7 + aparam(30)*x6y + &
          aparam(31)*x5y2 + aparam(32)*x4y3 + aparam(33)*x3y4 + &
          aparam(34)*x2y5 + aparam(35)*xy6 + aparam(36)*y7
  end if
  if(mfit.ge.45) then
     x8 = x7*x
     x7y = x7*y
     x6y2 = x6*y2
     x5y3 = x5*y3
     x4y4 = x4*y4
     x3y5 = x3*y5
     x2y6 = x2*y6
     xy7 = x*y7
     y8 = y7*y
     surfaceval = surfaceval + aparam(37)*x8 + aparam(38)*x7y + &
          aparam(39)*x6y2 + aparam(40)*x5y3 + aparam(41)*x4y4 + &
          aparam(42)*x3y5 + aparam(43)*x2y6 + aparam(44)*xy7  + &
          aparam(45)*y8
  end if
  if(mfit.ge.55) then
     x9 = x8*x
     x8y = x8*y
     x7y2 = x7*y2
     x6y3 = x6*y3
     x5y4 = x5*y4
     x4y5 = x4*y5
     x3y6 = x3*y6
     x2y7 = x2*y7
     xy8 = x*y8
     y9 = y8*y
     surfaceval = surfaceval + aparam(46)*x9 + aparam(47)*x8y + &
          aparam(48)*x7y2 + aparam(49)*x6y3 + aparam(50)*x5y4 + &
          aparam(51)*x4y5 + aparam(52)*x3y6 + aparam(53)*x2y7 + &
          aparam(54)*xy8  + aparam(55)*y9
  end if
  if(mfit.eq.66) then
     x10 = x9*x
     x9y = x9*y
     x8y2 = x8*y2
     x7y3 = x7*y3
     x6y4 = x6*y4
     x5y5 = x5*y5
     x4y6 = x4*y6
     x3y7 = x3*y7
     x2y8 = x2*y8
     xy9 = x*y9
     y10 = y9*y
     surfaceval = surfaceval + aparam(56)*x10 + aparam(57)*x9y + &
          aparam(58)*x8y2 + aparam(59)*x7y3 + aparam(60)*x6y4 + &
          aparam(61)*x5y5 + aparam(62)*x4y6 + aparam(63)*x3y7 + &
          aparam(64)*x2y8 + aparam(65)*xy9 + aparam(66)*y10
  end if

  fparam(1) = 1.0
  if(mfit.ge.3) then
     fparam(2) = x
     fparam(3) = y
  end if
  if(mfit.ge.6) then
     fparam(4) = x2
     fparam(5) = x*y
     fparam(6) = y2
  end if
  if(mfit.ge.10) then
     fparam(7) = x3
     fparam(8) = x2y
     fparam(9) = xy2
     fparam(10) = y3
  end if
  if(mfit.ge.15) then
     fparam(11) = x4
     fparam(12) = x3y
     fparam(13) = x2y2
     fparam(14) = xy3
     fparam(15) = y4
  end if
  if(mfit.ge.21) then
     fparam(16) = x5
     fparam(17) = x4y
     fparam(18) = x3y2
     fparam(19) = x2y3
     fparam(20) = xy4
     fparam(21) = y5
  end if
  if(mfit.ge.28) then
     fparam(22) = x6
     fparam(23) = x5y
     fparam(24) = x4y2
     fparam(25) = x3y3
     fparam(26) = x2y4
     fparam(27) = xy5
     fparam(28) = y6
  end if
  if(mfit.ge.36) then
     fparam(29) = x7
     fparam(30) = x6y
     fparam(31) = x5y2
     fparam(32) = x4y3
     fparam(33) = x3y4
     fparam(34) = x2y5
     fparam(35) = xy6
     fparam(36) = y7
  end if
  if(mfit.ge.45) then
     fparam(37) = x8
     fparam(38) = x7y
     fparam(39) = x6y2
     fparam(40) = x5y3
     fparam(41) = x4y4
     fparam(42) = x3y5
     fparam(43) = x2y6
     fparam(44) = xy7
     fparam(45) = y8
  end if
  if(mfit.ge.55) then
     fparam(46) = x9
     fparam(47) = x8y
     fparam(48) = x7y2
     fparam(49) = x6y3
     fparam(50) = x5y4
     fparam(51) = x4y5
     fparam(52) = x3y6
     fparam(53) = x2y7
     fparam(54) = xy8
     fparam(55) = y9
  end if
  if(mfit.eq.66) then
     fparam(56) = x10
     fparam(57) = x9y
     fparam(58) = x8y2
     fparam(59) = x7y3
     fparam(60) = x6y4
     fparam(61) = x5y5
     fparam(62) = x4y6
     fparam(63) = x3y7
     fparam(64) = x2y8
     fparam(65) = xy9
     fparam(66) = y10
  end if

end subroutine surface

!-----------------------------------------------------------------------

! Subroutine surface1(xx,aparam,mfit,dum1,dum2,> fparam,surface1val)
! Subroutine to do a quadratic fit to the magnitude differences.

subroutine surface1(xx,aparam,mfit,dum1,dum2,fparam,surface1val)

  real :: xx(2)  !2 dim only neccessary for coherence in subrout. chisq
  real :: aparam(mfit),fparam(mfit)
  real :: dum1,dum2    !Only neccessary for coherence in subrout. chisq

  x=xx(1)

  surface1val = aparam(mfit)
  fparam(1) = 1
  i = mfit-1
  do while (i.ge.1)
     surface1val = surface1val*x+aparam(i)
     fparam(mfit+1-i) = fparam(mfit-i)*x
     i = i-1
  end do

end subroutine surface1

!---------------------------------------------------------------------------
!
! Subroutine surfpoly(z,x,y,n,a,mfit)
!
! This subroutine does a fit to a surface using a linear least-squares routine.
! The linear least-squares fitting is done by solving the linear equations
! by Gauss-Jordan elimination

subroutine surfpoly(z,x,y,n,a,mfit)

  implicit real(8) (a-h,o-z)
  integer,parameter :: mmax=55
  real :: a(mmax)
  real :: z(n),x(n),y(n)
  real(8) :: alpha(mmax,mmax),beta(mmax)

!  We need to solve the matrix equation beta = a * alpha
!  First, build the matrices.
  alpha = 0
  beta = 0
  do i=1,n
     xx = dble(x(i))
     yy = dble(y(i))
     zz = dble(z(i))
     do j=1,mfit
        do k=1,mfit
           alpha(j,k) = alpha(j,k)+term(j,xx,yy)*term(k,xx,yy)
        end do
        beta(j) = beta(j)+term(j,xx,yy)*zz
     end do
  end do
! Solve this system using the GJ elimination routine (set to double precision)
  call gaussj(alpha,mfit,mmax,beta,1,1)

  do i=1,mfit
     a(i) = sngl(beta(i))      !Beta is transformed in gaussj into a
  end do

end subroutine surfpoly

!---------------------------------------------------------------------------
!
! Function surfeval(x,y,a,mfit)
! Function to evaluate the surface function solved above.

function surfeval(x,y,a,mfit)

  real(8) :: xx,yy,aa,sum,term
  real :: a(mfit)

  xx = dble(x)
  yy = dble(y)

  sum = 0.0
  do i=1,mfit
     aa = dble(a(i))
     sum = sum+aa*term(i,xx,yy)
  end do

  surfeval = sngl(sum)
  
end function surfeval

!---------------------------------------------------------------------------
!
! Function term(i,x,y)
! Function to evaluate the terms of the linear equations used above.

function term(i,x,y)

  implicit real(8) (a-h,o-z)

  if(i.gt.55) then
     term = 0.0
     return
  end if

  if(i.eq.1) term = 1.0

  if(i.eq.2) term = x
  if(i.eq.3) term = y

  if(i.eq.4) term = x**2
  if(i.eq.5) term = x*y
  if(i.eq.6) term = y**2

  if(i.eq.7) term = x**3
  if(i.eq.8) term = (x**2)*y
  if(i.eq.9) term = x*(y**2)
  if(i.eq.10) term = y**3

  if(i.eq.11) term = x**4
  if(i.eq.12) term = (x**3)*y
  if(i.eq.13) term = (x**2)*(y**2)
  if(i.eq.14) term = x*(y**3)
  if(i.eq.15) term = y**4

  if(i.eq.16) term = x**5
  if(i.eq.17) term = (x**4)*y
  if(i.eq.18) term = (x**3)*(y**2)
  if(i.eq.19) term = (x**2)*(y**3)
  if(i.eq.20) term = x*(y**4)
  if(i.eq.21) term = y**5

  if(i.eq.22) term = x**6
  if(i.eq.23) term = (x**5)*y
  if(i.eq.24) term = (x**4)*(y**2)
  if(i.eq.25) term = (x**3)*(y**3)
  if(i.eq.26) term = (x**2)*(y**4)
  if(i.eq.27) term = x*(y**5)
  if(i.eq.28) term = y**6

  if(i.eq.29) term = x**7
  if(i.eq.30) term = (x**6)*y
  if(i.eq.31) term = (x**5)*(y**2)
  if(i.eq.32) term = (x**4)*(y**3)
  if(i.eq.33) term = (x**3)*(y**4)
  if(i.eq.34) term = (x**2)*(y**5)
  if(i.eq.35) term = x*(y**6)
  if(i.eq.36) term = y**7

  if(i.eq.37) term = x**8
  if(i.eq.38) term = (x**7)*y
  if(i.eq.39) term = (x**6)*(y**2)
  if(i.eq.40) term = (x**5)*(y**3)
  if(i.eq.41) term = (x**4)*(y**4)
  if(i.eq.42) term = (x**3)*(y**5)
  if(i.eq.43) term = (x**2)*(y**6)
  if(i.eq.44) term = x*(y**7)
  if(i.eq.45) term = y**8

  if(i.eq.46) term = x**9
  if(i.eq.47) term = (x**8)*y
  if(i.eq.48) term = (x**7)*(y**2)
  if(i.eq.49) term = (x**6)*(y**3)
  if(i.eq.50) term = (x**5)*(y**4)
  if(i.eq.51) term = (x**4)*(y**5)
  if(i.eq.52) term = (x**3)*(y**6)
  if(i.eq.53) term = (x**2)*(y**7)
  if(i.eq.54) term = x*(y**8)
  if(i.eq.55) term = y**9

end function term

!---------------------------------------------------------------------------
!
! Subroutine gaussj(a,n,np,b,m,mp)
! Linear equation solution by Gauss-Jordan elimination.
! Notice that on the output, a is replaced by its matrix inverse, 
! and b is replaced by the corresponding set of solution vectors.
! To better understand this subroutine, check section 2.1 
! in Numerical Recipes in Fortran.
! Notice that this is a real(8) version of Num Rec gaussj.

subroutine gaussj(a,n,np,b,m,mp)

  implicit real(8) (a-h,o-z)
  integer, parameter :: nmax=55
  integer :: m,mp,n,np
  real(8) ::  a(np,np),b(np,mp)
  integer :: indxc(nmax),indxr(nmax),ipiv(nmax)
 
  do j=1,n
     ipiv(j)=0
  end do
  do i=1,n
     big=0.
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if (ipiv(k).eq.0) then
                 if (dabs(a(j,k)).ge.big)then
                    big=dabs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 print*,'singular matrix in gaussj'
                 stop
              endif
           enddo
        endif
     enddo
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) then
        print*,'singular matrix in gaussj'
        stop
     end if
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo
     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo
        endif
     enddo
  enddo
  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  enddo

end subroutine gaussj

!---------------------------------------------------------------------------
!
! Function mgetmfit(norder)
! This function gives the number of terms in a bivariate polynomial of n order

function mgetmfit(norder)

  integer :: mgetmfit

  mgetmfit = 0
  do i=1,norder+1
     mgetmfit = mgetmfit + i
  end do

end function mgetmfit


